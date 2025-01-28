
'''

NOTE: this code doesn't currently run. It relies on an old version of the code, which I didn't bother porting into the new version for running these experiments.
Functions like `load_hcnn_models` and `make_TCR_pMHC_predictions_from_pose` all have their respective equivalents in the new codebase.
I will update this as soon as possible.

'''

from argparse import ArgumentParser
import json
import logging
import os
from pathlib import Path
import sys
import time
from typing import Dict, List, Optional

from copy import deepcopy

from functools import partial
import h5py
import numpy as np
import pandas as pd
from multiprocessing import Pool
import pyrosetta
from scipy.special import softmax
from tqdm import tqdm
from typing import List, Tuple

from runtime.tcr_pmhc_experiments.peptide_generation_with_sim_anneal.temperature.multiplicative import (
    MultiplicativeScheduler)

from runtime.tcr_pmhc_experiments.peptide_generation_with_sim_anneal.temperature.thermodynamic import(
    ThermodynamicScheduler)

from runtime.tcr_pmhc_experiments.utils.predict import (
    load_hcnn_models
)

from protein_holography_pytorch.utils.hcnn.prediction_parser import (
    get_energy)

from protein_holography_pytorch.utils.protein_naming import (
    ind_to_ol_size, ol_to_ind_size)


from runtime.tcr_pmhc_experiments.utils.predict import (
    make_TCR_pMHC_predictions_from_pose)

from runtime.tcr_pmhc_experiments.utils.pocket import (
    get_tcr_pmhc_peptide_pocket_resnums, get_pockets)

from runtime.tcr_pmhc_experiments.utils.pep import (
    pdb_to_pep)


import protein_holography_pytorch.utils.pyrosetta_utils as pru

from protein_holography_pytorch.utils.pyrosetta_utils import (
    default_flags, wet_flags, mute_flag,
    get_pose_residue_number, get_pdb_residue_info,
    find_calpha_neighbors,
    fastrelax_positions, make_mutations)

# # for testing purposes
# from protein_holography.tcr_pmhc import pocket
# from runtime.tcr_pmhc_experiments.peptide_generation_with_sim_anneal.temperature import thermodynamic, multiplicative
# from importlib import reload
# reload(pocket)
# reload(thermodynamic)
# reload(multiplicative)
# get_pockets = pocket.get_pockets
# ThermoynamicsSchduler = thermodynamic.ThermodynamicScheduler
# MultiplicativeScheduler = multiplicative.MultiplicativeScheduler


N_AA = 20 # number of amino-acid types


class TCR_pMHC_Annealer():

    def __init__(
            self,
            mutant_tup: Tuple,
            peptide_L: int,
            T_regions: List,
            iterations: int,
            hcnn_config_file: str,
            region_to_optimize_energy_of: str,
            max_atoms: int=12000,
            peptide_backbone_relax: bool=True,
            peptide_acceptance: bool=False,
            force_best_mutations: bool=False,
            write_pdbdir: str=None,
            track_mutate_metrics: bool=False,
            testing_mode: bool=False,
            mhc_muts: Optional[List]=None,
            mhc_chain: str = 'A',
            use_logp_as_energy: bool = False
    ):
        """
        Parameters
        ----------
        
        """
        self.T_regions = T_regions 
        self.peptide_L = peptide_L
        self.max_atoms = max_atoms
        self.iterations = iterations
        self.load_network(hcnn_config_file)
        self.init_score_fxns()
        self.get_pdb(mutant_tup)
        self.track_mutate_metrics = track_mutate_metrics
        self.acceptance = peptide_acceptance
        self.force = force_best_mutations
        self.testing_mode = testing_mode
        self.peptide_backbone_relax = peptide_backbone_relax
        self.mhc_muts = mhc_muts
        self.mhc_chain = mhc_chain
        self.use_logp_as_energy = use_logp_as_energy
        self.energy_str = 'pnlogp' if self.use_logp_as_energy else 'pnE'
        self.region_to_optimize_energy_of = region_to_optimize_energy_of
        # establish schedulers for each region
        #
        # param dictionary should include the following pieces of info
        # in this order

        # T_0, T_f,
        self.establish_pose()

        logging.info(f"MHC mutations are {mhc_muts}")
        if mhc_muts != None:
            logging.info(f"Making mhc mutations {''.join(self.mhc_muts)}")
            self.make_mhc_mutations(self.mhc_muts, chain=self.mhc_chain)

        self.establish_history()

        self.schedulers = {}
        for region, params in T_regions.items():
            self.schedulers[region] = self.create_scheduler(**params)
            if self.schedulers[region] == None:
                logging.warning(f"Region {region} supplied for temperature "
                                "but schedule could not be constructed")

        if track_mutate_metrics:
            self.add_mutate_regions()

        self.write_pdb = (write_pdbdir != None)
        if self.write_pdb:
            self.write_pdbdir = write_pdbdir
        

    def make_mhc_mutations(self, mhc_muts: List[str], chain: str='A'):
        """Make mutations to MHC"""
        N_muts = len(mhc_muts)
        wt_aas, sites, mut_aas = [], [], []
        icodes = ['']*N_muts
        mhc_seqposs = []
        for mhc_mut in mhc_muts:
            wt_aa = mhc_mut[0]
            site = int(mhc_mut[1:-1])
            mut_aa = mhc_mut[-1]
            wt_aas.append(wt_aa)
            mut_aas.append(mut_aa)
            sites.append(site)
            mhc_seqpos = get_pose_residue_number(self.pose, chain, site)
            mhc_seqposs.append(mhc_seqpos)
        self.mutate_sites(mhc_seqposs, mut_aas)


    def load_network(self, hcnn_config_file):
        # load HCNN model, or suite of HCNN models
        with open(hcnn_config_file, 'r') as f:
            model_config = json.load(f)
        model_name = model_config['name']
        logging.info(f"Loading {model_name}")
        model_dir_list = model_config['model_dir_list']
        hcnn_models, hparams = load_hcnn_models(model_dir_list)
        self.hcnn_models = hcnn_models
        self.hparams = hparams

    def init_score_fxns(self):
        # scorefunctions
        self.scorefxn = pyrosetta.create_score_function('ref2015.wts')
        self.scorefxn_cartesian = pyrosetta.create_score_function('ref2015_cart.wts')

    def get_pdb(self, mutant_tup: Tuple):
        # extract info from mutant tup
        self.pdbfile = mutant_tup[0]
        self.pdb = Path(self.pdbfile).stem
        self.chain = mutant_tup[1]
        self.mutation_inds = mutant_tup[2]    
        self.init_seq = ''.join(mutant_tup[3])
        self.mutant_tup = mutant_tup

    def establish_pose(self):
        # make pose and get relaxation sites
        self.pose = pyrosetta.pose_from_pdb(self.pdbfile)
        
        self.peptide_resnums, self.pocket_resnums, self.complex_resnums = get_tcr_pmhc_peptide_pocket_resnums(self.pose, self.mutation_inds, self.chain)
        self.pocket_size = len(self.pocket_resnums)

        self.regions = get_pockets(self.pdb, self.pocket_resnums, self.peptide_resnums, self.pose, mhc_chain=self.mhc_chain)
                
    def create_scheduler(self, **params):
        print(params)
        scheduler = None
        if params["schedule"] == "multiplicative":
            logging.info("Creating multiplicative scheduler")
            scheduler = MultiplicativeScheduler(
                metric_history=self.metric_history, **params
            )
        if params["schedule"] == "thermodynamic":
            logging.info("Creating thermodynamic scheduler")
            scheduler = ThermodynamicScheduler(
                metric_history=self.metric_history,
                **params)

        return scheduler

    def add_mutate_regions(self):
        self.mutate_regions = {}
        for region,resnums in self.regions.items():
            self.mutate_regions[f'mutate_{region}'] = resnums

    def establish_history(self):
        """
        Create metric_history dictionary to keep track of different 
        observables throughout the annealing process
        """
        self.metric_history = {}
        sequences = np.zeros(shape=(self.iterations + 1,), dtype='S11')
        self.metric_history["sequences"] = sequences
        for region, seqnums in self.regions.items():
            region_length = len(seqnums)
            self.metric_history.update(
                {
                    f'{region} pes': np.zeros(
                        shape=(self.iterations, region_length, N_AA), dtype=float),
                    f'{region} pnE': np.zeros(
                        shape=(self.iterations,), dtype=float),
                    f'{region} pnlogp': np.zeros(
                        shape=(self.iterations,), dtype=float)
                }
            )
        if self.track_mutate_metrics:
            for region, seqnums in self.mutate_regions.items():
                region_length = len(seqnums)
                self.metric_history.update(
                    {
                    f'{region} pes': np.zeros(
                        shape=(self.iterations, region_length, N_AA), dtype=float),
                    f'{region} pnE': np.zeros(
                        shape=(self.iterations,), dtype=float),
                    f'{region} pnlogp': np.zeros(
                        shape=(self.iterations,), dtype=float)
                    }
                )
        self.metric_history['acceptance'] = np.zeros(
            shape=(self.iterations,), dtype=int)
        self.metric_history[f"next_peptide {self.energy_str}"] = np.zeros(
            shape=(self.iterations,), dtype=float)
        for region in self.T_regions.keys():
            self.metric_history[f"{region} temperature"] = np.zeros(
                shape=(self.iterations,), dtype=float)


    def mutate_sites(self, sites: List[int], mut_aas: List[str]):
        mutations = {
            seqpos:mutation for seqpos,mutation in
            zip(sites, mut_aas)
        }
        logging.info(f"Making initial mutations to {mut_aas}")
        make_mutations(mutations, self.pose)

    # def mutate(self):
    #     mutations = {
    #         seqpos:mutation for seqpos,mutation in
    #         zip(self.peptide_resnums, self.mutant_tup[3])
    #     }
    #     logging.info(f"Making initial mutations to {self.mutant_tup[3]}")
    #     make_mutations(mutations, self.pose)
            
    # def relax(self):
    #     start_time = time.time()
    #     if self.testing_mode:
    #         logging.info("relaxing structure in testing mode")
    #         fastrelax_positions(
    #             self.scorefxn_cartesian,
    #             [], 
    #             self.peptide_resnums, 
    #             self.pose
    #         )
    #     else:
    #         logging.info("relaxing structure")
    #         fastrelax_positions(
    #             self.scorefxn_cartesian,
    #             self.peptide_resnums,
    #             self.complex_resnums, #pocket_resnums,
    #             self.pose
    #         )
    #     logging.debug("--- %s second for fastrelax---" % (time.time() - start_time))

    def get_original_pocket_energy(self):
        results = make_TCR_pMHC_predictions_from_pose(
            self.pose,
            self.hcnn_models,
            self.hparams,
            self.regions,
            average_metrics=True,
            compute_zgrams_only_for_requested_regions=True
        )
        return results

    def check_pocket_energies(
        self,
        beta: float,
        i: int
    ) -> bool:

        accept = True
        reject = False
        
        # check pocket energy
        new_pocket_energy = self.metric_history[f'{self.region_to_optimize_energy_of} {self.energy_str}'][i]
        if i == 0:
            old_pocket_energy = 0.
        else:
            old_pocket_energy = self.metric_history[f'{self.region_to_optimize_energy_of} {self.energy_str}'][self.metric_history['acceptance'] == 1][-1]
        logging.info(f"Old {self.region_to_optimize_energy_of} energy {old_pocket_energy}")
        logging.info(f"New {self.region_to_optimize_energy_of} energy {new_pocket_energy}")
        logging.info(f"{self.region_to_optimize_energy_of} energies {self.metric_history[f'{self.region_to_optimize_energy_of} {self.energy_str}']}")
        if new_pocket_energy < old_pocket_energy:
            delta_E = old_pocket_energy - new_pocket_energy
            logging.debug(f"delta E = {delta_E}")
            if np.random.uniform(size=(1,))[0] < np.exp(-1 * beta * (old_pocket_energy - new_pocket_energy)):
                logging.info(f"New {self.region_to_optimize_energy_of} is worse but sampling anyway")
                return accept
            else:
                logging.info(f"New {self.region_to_optimize_energy_of} is worse and rejecting")
                return reject
        logging.info(f"New {self.region_to_optimize_energy_of} is better. Accepting automatically")
        return accept
                
    def anneal(self):
        """Anneal the structure"""
        logging.info("Beginning annealing process")
        self.metric_history["sequences"][0] = self.init_seq
        mutations = {
            seqpos:mutation for seqpos,mutation in
            zip(self.peptide_resnums, self.mutant_tup[3])
        }
        logging.info(f"Making initial mutations to {self.mutant_tup[3]}")
        make_mutations(mutations, self.pose)
        k=0
        if self.write_pdb:
            os.chdir(self.write_pdbdir)
            print("dumping mutate cif")
            self.pose.dump_cif(
                os.path.join(f'step_{k}_mutate.cif')
            )

        logging.info("relaxing structure")
        logging.debug(f"Backbone relax sites are {self.peptide_resnums * self.peptide_backbone_relax}")
        fastrelax_positions(
            self.scorefxn_cartesian,
            self.peptide_resnums * self.peptide_backbone_relax,
            self.complex_resnums, # <-- complex_resnums is the full pocket within 10. AA of the peptide
            self.pose
        )
        if self.write_pdb:
            os.chdir(self.write_pdbdir)
            self.pose.dump_cif(
                os.path.join(f'step_{k}_relax.cif')
            )
        logging.info("Calculating energetics")
        results = make_TCR_pMHC_predictions_from_pose(
            self.pose, self.hcnn_models, self.hparams, self.regions,
            #resnums=(peptide_resnums,pocket_resnums), # <- to change
            average_metrics=True,
            compute_zgrams_only_for_requested_regions=True
        )

        last_pose = self.pose.clone() # NOTE: this might not be doing anything!!! just a pass-by-reference that will change with self.pose
        last_results = deepcopy(results)
        last_mutant_tup = deepcopy(self.mutant_tup)

        # establish initial temperatures
        peptide_beta = 1./self.schedulers["peptide"].T_0
        pocket_beta = 1./self.schedulers["pocket"].T_0

        self.metric_history["sequences"][0] = ''.join(self.mutant_tup[3])
        self.metric_history["acceptance"][0] = 1
        print(self.metric_history.keys())
        print(results.keys())
        for region in self.regions.keys():
            self.metric_history[f'{region} pes'][0] = results[f'{region} pes']
            self.metric_history[f'{region} pnE'][0] = results[f'{region} pnE']
            self.metric_history[f'{region} pnlogp'][0] = results[f'{region} pnlogp']
        try:
            self.metric_history[f'complex pnE'][0] = results[f'pocket pnE'] + results[f'peptide pnE']
            self.metric_history[f'complex pnlogp'][0] = results[f'pocket pnlogp'] + results[f'peptide pnlogp']
        except Exception as e:
            print(e)
            print('Available regions: ', results.keys())
            pass
        
        for k in tqdm(range(1, self.iterations)):

            print('-'*40)
            print(f'Address of last_pose:\t{id(last_pose)}')
            print(f'Address of self.pose:\t{id(self.pose)}')
            print('-'*40)

            logging.info(f"Sampling mutant at step {k}")

            mutations_ol = sample_mutant(results["peptide pes"],
                                         peptide_beta,
                                         acceptance=self.acceptance,
                                         force=self.force)
            
            logging.info(f"Next sequence is {''.join(mutations_ol)}")

            new_network_energy = get_energy(''.join(mutations_ol), results["peptide pes"])
            self.metric_history[f"next_peptide {self.energy_str}"][k-1] = new_network_energy
            
            logging.debug(f'new network energy: {new_network_energy}')
            logging.debug(f'curr struct network energy: {results[f"peptide {self.energy_str}"]}')


            # # check peptide acceptance condition
            # if self.acceptance:
            #     logging.info('Checking acceptance crietrion')
            #     # if neq seq is better, accept by default
            #     if new_network_energy > curr_network_energy: 
            #         mutant_tup = (*mutant_tup[:-1],mutations_ol)
            #         curr_network_energy = new_network_energy
            #         # if new seq worse, accept with some probability
            #     elif np.random.uniform(
            #             size=(1,))[0] < softmax(
            #                 beta * np.array([
            #                     curr_network_energy,new_network_energy]))[1]:
            #         mutant_tup = (*mutant_tup[:-1],mutations_ol)
            #         curr_network_energy = new_network_energy
            #         logging.info('New sequence is worse but sampling anyway')
            #     else: # new sequence rejected
            #         logging.info('Rejecting new sequence')
            # else: # no acceptance criteria. Automatically accept new sequence
            #     mutant_tup = (*mutant_tup[:-1],mutations_ol)

            self.mutant_tup = (*self.mutant_tup[:-1], mutations_ol)
            
            logging.info(f"Making mutations to {mutations_ol}")
            mutations = {
                seqpos:mutation for seqpos,mutation in
                zip(self.peptide_resnums, self.mutant_tup[3])
            }

            # print('-'*80)
            # print(f'mutations_ol:\t{mutations_ol}')
            
            make_mutations(mutations, self.pose)

            if self.write_pdb:
                os.chdir(self.write_pdbdir)
                print("dumping mutate cif")
                self.pose.dump_cif(
                    os.path.join(f'step_{k}_mutate.cif')
                )

            logging.info("Relaxing structure")
            fastrelax_positions(
                self.scorefxn_cartesian,
                self.peptide_resnums * self.peptide_backbone_relax,
                self.complex_resnums,  # <-- complex_resnums is the full pocket within 10. AA of the peptide
                self.pose
            )

            if self.write_pdb:
                os.chdir(self.write_pdbdir)
                self.pose.dump_cif(
                    os.path.join(f'step_{k}_relax.cif')
                )
            
            results = make_TCR_pMHC_predictions_from_pose(
                self.pose, self.hcnn_models, self.hparams, self.regions,
                #resnums=(peptide_resnums,pocket_resnums), # <- to change
                average_metrics=True,
                compute_zgrams_only_for_requested_regions=True
            )

            self.metric_history["sequences"][k] = ''.join(self.mutant_tup[3]) 
            for region in self.regions.keys():
                self.metric_history[f'{region} pes'][k] = results[f'{region} pes']
                self.metric_history[f'{region} pnE'][k] = results[f'{region} pnE']
                self.metric_history[f'{region} pnlogp'][k] = results[f'{region} pnlogp']
            try:
                self.metric_history[f'complex pnE'][k] = results[f'pocket pnE'] + results[f'peptide pnE']
                self.metric_history[f'complex pnlogp'][k] = results[f'pocket pnlogp'] + results[f'peptide pnlogp']
            except Exception as e:
                print(e)
                print('Available regions: ', results.keys())
                pass
                
            logging.info("Checking pocket acceptance")
            accept = self.check_pocket_energies(pocket_beta, k)

            if not accept:
                self.pose = last_pose
                results = last_results # NOTE - FIX: this changed self.results -> results, which makes things different
                self.mutant_tup = last_mutant_tup # this doesn't matter really since it gets resampled everytiome from "results"
            else:
                self.metric_history["acceptance"][k] = 1
                ## NOTE - FIX: these 4 lines have been added, and might change things
                # from copy import deepcopy
                last_pose = self.pose.clone() # NOTE: this isn't actually a full deep copy, would need to use detached_copy() but need to learn the syntax (takes as input a pose object )
                last_results = deepcopy(results)
                last_mutant_tup = deepcopy(self.mutant_tup) # this doesn't matter really since it gets resampled everytiome from "results"

            logging.info("Stepping temperature schedulers")
            # step schedulers and calculate next temperatures
            for scheduler in self.schedulers.values():
                scheduler.step()
            logging.info("pep temps",self.metric_history["peptide temperature"])
            logging.info("pocket temps", self.metric_history["pocket temperature"])
            if k != self.iterations - 1:
                peptide_beta = 1./self.metric_history["peptide temperature"][k+1]
                pocket_beta = 1./self.metric_history["pocket temperature"][k+1]
                logging.info(f"Next peptide T = {1./peptide_beta}")
                logging.info(f"Next pocket T = {1./pocket_beta}")
        
        logging.info("All annealing steps run successsfully")


def sample_mutant(
    pseudoenergies: np.ndarray, 
    beta: float,
    acceptance: bool=False,
    force: bool=False
) -> List:
    probs = softmax(beta * pseudoenergies.astype(np.float64),axis=-1)

    mutations = np.zeros(shape=probs.shape)
    for j,prob in enumerate(probs):
        if force:
            logging.debug('Forcing most likely mutants')
            curr_state = ol_to_ind_size[mutant_tup[-1][j]]
            if curr_state != np.argmax(prob):
                mutations[j] = np.zeros(shape=(20))
                mutations[j,np.argmax(prob)] = 1
                state_ol = ind_to_ol_size[np.argmax(prob)]
                logging.info(
                    f'Forcing acceptance of state '
                    f'{state_ol} at site {j}')
            else:
                try:
                    mutations[j] = np.random.multinomial(1,prob/np.sum(prob))
                except Exception as e:
                    print('Error with multinomial sampling:', e)
                    print(prob)
                    print(np.sum(prob[:-1]))
                    print(np.sum(prob))
                    print(np.sum(prob[:-1]/np.sum(prob)))
                    print(np.sum(prob/np.sum(prob)))
                    mutations[j] = np.zeros(shape=(20))
                    mutations[j,np.argmax(prob)] = 1
        else:
            try:
                mutations[j] = np.random.multinomial(1, prob/np.sum(prob))
            except Exception as e:
                print('Error with multinomial sampling:', e)
                print(prob)
                print(np.sum(prob[:-1]))
                print(np.sum(prob))
                print(np.sum(prob[:-1]/np.sum(prob)))
                print(np.sum(prob/np.sum(prob)))
                mutations[j] = np.zeros(shape=(20))
                mutations[j,np.argmax(prob)] = 1
                
    mutations_ol = [ind_to_ol_size[x] for x in np.argmax(mutations,axis=-1)]

    return mutations_ol


if __name__=="__main__":
    parser = ArgumentParser()
    parser.add_argument('--pdb', type=str, required=True)
    parser.add_argument('--pdbdir',type=str, default='/gscratch/spe/gvisan01/tcr_pmhc/pdbs')
    parser.add_argument('--output_dir', type=str, default='/gscratch/scrubbed/gvisan01/tcrpmhc/annealing_runs')

    parser.add_argument('--chain', type=str, default='C',
                        help='Peptide chain to mutate')
    
    parser.add_argument('--peptide_resnum_start', type=int, default=1)
    
    parser.add_argument('--mhc_chain', type=str, default='A',
                        help='MHC chain, for the purpose of mutating it (if inputting mhc_muts) and keeping track of its specific energy.')
    
    parser.add_argument('--mhc_muts', default=None, type=str, nargs='+',
                        help='MHC mutations to make, in standard [WT][site][mut] format. e.g. A9V.')
    
    parser.add_argument('--sequence', type=str, default='',
                        help='Initial peptide sequence for annealing. Empty string for random sequence')
        
    parser.add_argument('--schedule_file', type=str, required=True,
                        help='Path to JSON file containing annealing schedule')

    parser.add_argument('--iters', type=int, default=100)

    parser.add_argument('--energy_to_use', type=str, default='pnE', choices=['pnE', 'pnlogp', 'dsmbind_pp'])
    parser.add_argument('--region_to_optimize_energy_of', type=str, default='pocket', choices=['peptide', 'pocket', 'complex'])

    parser.add_argument('--max_atoms', type=int, default=15000)

    parser.add_argument('--peptide_acceptance', action='store_true', default=False)
    parser.add_argument('--force_best_mutations', action='store_true', default=False)
    parser.add_argument('--add_crystal_water', action="store_true", default=False)

    parser.add_argument('--dont_relax_peptide_backbone', dest='relax_peptide_backbone', action='store_false', default=True)

    parser.add_argument('--model_config', type=str, required=True,
                        help='Path to json file that contains relevant config values for the model(s) to use. For now, contains an identifier and a list of model_dirs.')

    parser.add_argument('--write_pdb', action='store_true', default=False)

    parser.add_argument('--logging', default='INFO', type=str)
    parser.add_argument('--pyrosetta_logging', default='WARNING')
    parser.add_argument('--rosetta_logging', default='WARNING')

    parser.add_argument('--job', type=int, default=0,
                        help='Dummy job number for ease of slurm job submission')

    args = parser.parse_args()


    logging.getLogger().setLevel(args.logging)
    if args.pyrosetta_logging != None:
        logging.getLogger("pyrosetta").setLevel(args.pyrosetta_logging)
    if args.rosetta_logging != None:
        logging.getLogger("rosetta").setLevel(args.rosetta_logging)

    logging.info(f"Reading schedule file {args.schedule_file}")
    with open(f'./temperature/schedule_files/{args.schedule_file}.json', 'r') as f:
        schedules = json.load(f)
        
    # prepare protein pose    
    pdbpath = os.path.join(args.pdbdir, args.pdb+'.pdb')
    if args.add_crystal_water:
        pyrosetta.init(wet_flags + mute_flag, silent=True)
        pose = pyrosetta.pose_from_pdb(pdbpath)
        wet_str = 'wet'
        logging.info('Wet protein made')
    else:
        pyrosetta.init(default_flags + mute_flag, silent=True)
        pose = pyrosetta.pose_from_pdb(pdbpath)
        wet_str = 'dry'
        logging.info('Dry protein made')
    
    # load HCNN model, or suite of HCNN models
    with open(args.model_config, 'r') as f:
        model_config = json.load(f)
    model_name = model_config['name']
        

    # get peptide residue ids

    N_mutants = 1
    if args.sequence != '':
        peptide_L = len(args.sequence)
    else:
        peptide_L = len(pdb_to_pep[args.pdb])

    peptide_seqnums = []
    for i in range(args.peptide_resnum_start, peptide_L + args.peptide_resnum_start):
        seqnum = pru.get_pose_residue_number(pose,args.chain,i)
        if seqnum != 0:
            peptide_seqnums.append(seqnum)

    ids = [pru.get_pdb_residue_info(pose,i) for i in peptide_seqnums]
    
    if args.sequence == "":
        sequences_int = np.random.randint(low=0,high=19,size=(N_mutants*peptide_L))
        sequences_ol = [ind_to_ol_size[x] for x in sequences_int]
    else:
        sequences_ol = [x for x in args.sequence]
    sequences_ol = np.reshape(sequences_ol,(N_mutants,peptide_L))
        
    # actual assignment of mutants
    mutant_tups = [
        (pdbpath,
         [x[0] for x in ids],[x[1] for x in ids],seq)
         for i,(temp_pdb,seq) in enumerate(zip([pdbpath],sequences_ol))
    ]
    mutant_tup = mutant_tups[0]
    original_mutant = ''.join(mutant_tup[-1])
    logging.info(f"original mutant is {original_mutant}")

    
    logging.info('Making pdb dir')
    schedule_name = Path(args.schedule_file).stem
    pdb_dir = f'anneal__{args.pdb}__{original_mutant}__{model_name}__{args.iters}__{schedule_name}__{args.energy_to_use}__{args.region_to_optimize_energy_of}__{wet_str}__{args.job}__pdbs'
    history_dir = pdb_dir.replace('__pdbs', '')
    os.makedirs(args.output_dir, exist_ok=True)
    
    logging.info(f'acceptance status = {args.peptide_acceptance}')
    logging.info(f'force status = {args.force_best_mutations}')

    annealer = TCR_pMHC_Annealer(
        mutant_tup=mutant_tups[0],
        peptide_L=peptide_L,
        T_regions=schedules,
        iterations=args.iters,
        hcnn_config_file=args.model_config,
        region_to_optimize_energy_of=args.region_to_optimize_energy_of,
        peptide_backbone_relax=args.relax_peptide_backbone,
        peptide_acceptance=args.peptide_acceptance,
        force_best_mutations=args.force_best_mutations, 
        write_pdbdir=(os.path.join(args.output_dir, pdb_dir)) if args.write_pdb else None,
        mhc_muts=args.mhc_muts,
        use_logp_as_energy=(args.energy_to_use=='pnlogp')
    )

    # just if we want to get the energy of the wildtype in the pdb structure
    results = annealer.get_original_pocket_energy()
    np.save(os.path.join(args.output_dir, f'original_pocket_energies__{args.pdb}__{model_name}__{wet_str}.npy'), results)

    from time import time
    start = time()
    annealer.anneal()
    end = time()
    print(f"Annealing took {end-start} seconds")
    
    logging.info('Saving history')
    np.save(os.path.join(args.output_dir, history_dir + f'.npy'), annealer.metric_history)

