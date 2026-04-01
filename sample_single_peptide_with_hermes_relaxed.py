

'''
NOTE:
    1. Cannot have resnum holes nor icodes in the peptide chain
'''

from argparse import ArgumentParser
import json
import logging
import os
from pathlib import Path
import time
from copy import deepcopy

import numpy as np
import pyrosetta
from pyrosetta import Pose
from scipy.special import softmax, log_softmax
from tqdm import tqdm
from typing import List, Tuple

from temperature.scheduler import MultiplicativeScheduler, ThermodynamicScheduler

from hermes.utils.protein_naming import ind_to_ol_size, ol_to_ind_size

from hermes.pyrosetta_utils import (
    init_flags,
    find_calpha_neighbors,
    get_pose_residue_number,
    get_pdb_residue_info,
    fastrelax_positions,
    make_mutations
)

from hermes.inference.inference_hermes import predict_from_pdbfile, load_hermes_models


N_AA = 20 # number of amino-acid types

RADIUS = 10.0

def get_tcr_pmhc_peptide_pocket_resnums(
    pose: Pose,
    mutation_inds: List,
    chain: str    
) -> Tuple[List]:
    """Get TCR pMHC pocket resnums"""
    peptide_resnums = [
        get_pose_residue_number(pose, c, x) for x,c
        in zip(mutation_inds,chain)]
    complex_resnums = find_calpha_neighbors(peptide_resnums, RADIUS, pose)
    pocket_resnums = np.setdiff1d(
        list(complex_resnums), peptide_resnums)
    return peptide_resnums, pocket_resnums, complex_resnums

def get_energy(
    seq: str,
    energies: np.ndarray # pes or logps
) -> float:
    """Calculate pseudoenergy given a prediction matrix and a sequence"""

    seq_idxs = np.array([ol_to_ind_size[x] for x in seq])
    network_energy = np.sum(energies[np.arange(len(seq_idxs)), seq_idxs])

    return network_energy

def predict_from_pose_average_and_format_results(
    pose,
    hermes_models,
    hparams,
    batch_size,
    regions
):
    results_raw = predict_from_pdbfile(pose, hermes_models, hparams, batch_size, regions=regions, chain=None, add_same_noise_level_as_training=False, ensemble_with_noise=False)

    results = {}
    for region in results_raw.keys():
        sequence = ''.join(results_raw[region]['res_ids'].astype(str)[:, 0])

        pes = np.mean(results_raw[region]['logits'], axis=0)
        logps = log_softmax(pes.astype(np.float64), axis=1)

        pnE = get_energy(sequence, pes)
        pnlogp = get_energy(sequence, logps)

        results[f'{region} pes'] = pes
        results[f'{region} pnE'] = pnE
        results[f'{region} pnlogp'] = pnlogp

    return results


class TCR_pMHC_Annealer():

    def __init__(
            self,
            hermes_dir: str,
            model_version: str,
            mutant_tup: Tuple,
            peptide_L: int,
            T_regions: List,
            iterations: int,
            region_to_optimize_energy_of: str,
            max_atoms: int=12000,
            peptide_backbone_relax: bool=True,
            peptide_acceptance: bool=False,
            force_best_mutations: bool=False,
            write_pdbdir: str=None,
            track_mutate_metrics: bool=False,
            testing_mode: bool=False,
            use_logp_as_energy: bool = False,
            batch_size: int = 32,
    ):
        """
        Parameters
        ----------
        
        """
        self.hermes_dir = hermes_dir
        self.T_regions = T_regions 
        self.peptide_L = peptide_L
        self.max_atoms = max_atoms
        self.iterations = iterations
        self.model_version = model_version
        self.load_network()
        self.init_score_fxns()
        self.get_pdb(mutant_tup)
        self.track_mutate_metrics = track_mutate_metrics
        self.acceptance = peptide_acceptance
        self.force = force_best_mutations
        self.testing_mode = testing_mode
        self.peptide_backbone_relax = peptide_backbone_relax
        self.use_logp_as_energy = use_logp_as_energy
        self.energy_str = 'pnlogp' if self.use_logp_as_energy else 'pnE'
        self.region_to_optimize_energy_of = region_to_optimize_energy_of
        self.batch_size = batch_size
        # establish schedulers for each region
        #
        # param dictionary should include the following pieces of info
        # in this order

        # T_0, T_f,
        self.establish_pose()

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

    def load_network(self):
        trained_models_path = os.path.join(self.hermes_dir, 'trained_models', self.model_version)
        model_dir_list = [os.path.join(trained_models_path, model_rel_path) for model_rel_path in os.listdir(trained_models_path)]
        models, hparams, _ = load_hermes_models(model_dir_list)

        self.hermes_models = models
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

        self.regions = {
            'peptide': [get_pdb_residue_info(pose, i) for i in self.peptide_resnums],
            'pocket': [get_pdb_residue_info(pose, i) for i in self.pocket_resnums],
            'complex': [get_pdb_residue_info(pose, i) for i in self.complex_resnums]
        }

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
        make_mutations(self.pose, mutations)


    def get_original_pocket_energy(self):
        return predict_from_pose_average_and_format_results(self.pose, self.hermes_models, self.hparams, self.batch_size, self.regions)


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
        make_mutations(self.pose, mutations)
        k=0
        if self.write_pdb:
            logging.info("dumping mutate cif")
            self.pose.dump_cif(
                os.path.join(self.write_pdbdir, f'step_{k}_mutate.cif')
            )

        logging.info("relaxing structure")
        logging.debug(f"Backbone relax sites are {self.peptide_resnums * self.peptide_backbone_relax}")
        fastrelax_positions(
            self.pose,
            self.peptide_resnums * self.peptide_backbone_relax,
            self.complex_resnums, # <-- complex_resnums is the full pocket within 10 AA of the peptide, including the peptide
            self.scorefxn_cartesian,
        )

        if self.write_pdb:
            logging.info("dumping relax cif")
            self.pose.dump_cif(
                os.path.join(self.write_pdbdir, f'step_{k}_relax.cif')
            )
        logging.info("Calculating energetics")

        results = predict_from_pose_average_and_format_results(self.pose, self.hermes_models, self.hparams, self.batch_size, self.regions)

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
        
        for k in tqdm(range(1, self.iterations)):

            # print('-'*40)
            # print(f'Address of last_pose:\t{id(last_pose)}')
            # print(f'Address of self.pose:\t{id(self.pose)}')
            # print('-'*40)

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

            self.mutant_tup = (*self.mutant_tup[:-1], mutations_ol)
            
            logging.info(f"Making mutations to {mutations_ol}")
            mutations = {
                seqpos:mutation for seqpos,mutation in
                zip(self.peptide_resnums, self.mutant_tup[3])
            }
            
            make_mutations(self.pose, mutations)

            if self.write_pdb:
                logging.info("dumping mutate cif")
                self.pose.dump_cif(
                    os.path.join(self.write_pdbdir, f'step_{k}_mutate.cif')
                )

            logging.info("Relaxing structure")
            fastrelax_positions(
                self.pose,
                self.peptide_resnums * self.peptide_backbone_relax,
                self.complex_resnums, # <-- complex_resnums is the full pocket within 10 AA of the peptide, including the peptide
                self.scorefxn_cartesian,
            )

            if self.write_pdb:
                logging.info("dumping relax cif")
                self.pose.dump_cif(
                    os.path.join(self.write_pdbdir, f'step_{k}_relax.cif')
                )
            
            results = predict_from_pose_average_and_format_results(self.pose, self.hermes_models, self.hparams, self.batch_size, self.regions)


            self.metric_history["sequences"][k] = ''.join(self.mutant_tup[3]) 
            for region in self.regions.keys():
                self.metric_history[f'{region} pes'][k] = results[f'{region} pes']
                self.metric_history[f'{region} pnE'][k] = results[f'{region} pnE']
                self.metric_history[f'{region} pnlogp'][k] = results[f'{region} pnlogp']
                
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

    parser.add_argument('--hermes_dir', type=str, required=True,
                        help='Path to HERMES directory containing the `trained_models` directory with the model you want to use for inference.')
    
    parser.add_argument('-m', '--model_version', type=str, required=True,
                        help='Name of HERMES model you want to use. \
                              Use `hermes_py_000` or `hermes_bp_000` for the model trained without added noise to the atoms, \
                              `hermes_py_050` or `hermes_bp_050` for the model trained with 0.50 Angstrom noise to the atoms. \
                              Models with `_py_` in the name use pyrosetta to parse protein structures and were the ones used in the paper, \
                              whereas models with `_bp_` in the name use biopython.')

    parser.add_argument('-p', '--pdb_file', type=str, required=True,
                        help='Path to the PDB file containing the TCR-pMHC structure.')
    
    parser.add_argument('-o', '--output_file', type=str, required=True,
                        help='Path to the output .txt file where sampled peptide sequences will be saved (one sequence per line).')

    parser.add_argument('-c', '--peptide_chain', type=str, required=True,
                        help='Chain identifier for the peptide in the PDB file.')
            
    parser.add_argument('--initial_sequence', type=str, default='',
                        help='Initial peptide sequence for annealing. Empty string for random sequence')
        
    parser.add_argument('--schedule_name', type=str, required=True, 
                        help='Annealing schedule JSON file name (without the .json extension) as found in `./temperature/schedule_files`')

    parser.add_argument('--iters', type=int, default=100,
                        help='Number of simulated annealing iterations. Default: 100.\
                              NOTE: the schedules in this directory are designed for 100 iterations.')

    parser.add_argument('--energy_to_use', type=str, default='pnE', choices=['pnE', 'pnlogp'])

    parser.add_argument('--region_to_optimize_energy_of', type=str, default='peptide', choices=['peptide', 'pocket', 'complex'])

    parser.add_argument('--dont_relax_peptide_backbone', dest='relax_peptide_backbone', action='store_false', default=True)

    parser.add_argument('--save_metric_history', action='store_true', default=False,
                        help='Whether to save the metric history of the simulated annealing algorithm. It will be saved in a .npy file by the same name as `output_file`.')

    parser.add_argument('--write_pdb', action='store_true', default=False,
                        help='Whether to save the relaxed PDB files with each peptide. They will be saved in a directory identified by the same name as `output_file` and suffix `__pdbs`.')

    parser.add_argument('--max_atoms', type=int, default=15000)
    parser.add_argument('--peptide_acceptance', action='store_true', default=False)
    parser.add_argument('--force_best_mutations', action='store_true', default=False)
    parser.add_argument('--add_crystal_water', action="store_true", default=False)
    parser.add_argument('--logging', default='WARNING', type=str)
    parser.add_argument('--pyrosetta_logging', default='WARNING')
    parser.add_argument('--rosetta_logging', default='WARNING')

    args = parser.parse_args()


    logging.getLogger().setLevel(args.logging)
    if args.pyrosetta_logging != None:
        logging.getLogger("pyrosetta").setLevel(args.pyrosetta_logging)
    if args.rosetta_logging != None:
        logging.getLogger("rosetta").setLevel(args.rosetta_logging)

    logging.info(f"Reading schedule file {args.schedule_name}")
    with open(f'./temperature/schedule_files/{args.schedule_name}.json', 'r') as f:
        schedules = json.load(f)
    
    pdbid = Path(args.pdb_file).stem
        
    # prepare protein pose
    if args.add_crystal_water:
        pyrosetta.init(init_flags + ' -ignore_waters 0', silent=True)
        pose = pyrosetta.pose_from_pdb(args.pdb_file)
        wet_str = 'wet'
        logging.info('Wet protein made')
    else:
        pyrosetta.init(init_flags, silent=True)
        pose = pyrosetta.pose_from_pdb(args.pdb_file)
        wet_str = 'dry'
        logging.info('Dry protein made')
            

    # get peptide residue ids

    def chain_info(pose, chain_id):
        for i in range(1, pose.total_residue() + 1):
            if pose.pdb_info().chain(i) == chain_id:
                first_resnum = pose.pdb_info().number(i)
                break
        count = sum(1 for i in range(1, pose.total_residue() + 1)
                    if pose.pdb_info().chain(i) == chain_id)
        return count, first_resnum

    peptide_L, peptide_resnum_start = chain_info(pose, args.peptide_chain)

    N_mutants = 1
    if args.initial_sequence != '':
        if len(args.initial_sequence) != peptide_L:
            raise ValueError(f"Length of input sequence {len(args.initial_sequence)} does not match peptide length {peptide_L} in pdb structure.")

    peptide_seqnums = []
    for i in range(peptide_resnum_start, peptide_L + peptide_resnum_start):
        seqnum = get_pose_residue_number(pose,args.peptide_chain,i)
        if seqnum != 0:
            peptide_seqnums.append(seqnum)

    ids = [get_pdb_residue_info(pose,i) for i in peptide_seqnums]
    
    if args.initial_sequence == "":
        sequences_int = np.random.randint(low=0,high=19,size=(N_mutants*peptide_L))
        sequences_ol = [ind_to_ol_size[x] for x in sequences_int]
    else:
        sequences_ol = [x for x in args.initial_sequence]
    sequences_ol = np.reshape(sequences_ol,(N_mutants,peptide_L))
        
    # actual assignment of mutants
    mutant_tups = [
        (args.pdb_file,
         [x[0] for x in ids],[x[1] for x in ids],seq)
         for i,(temp_pdb,seq) in enumerate(zip([args.pdb_file],sequences_ol))
    ]
    mutant_tup = mutant_tups[0]
    original_mutant = ''.join(mutant_tup[-1])
    logging.info(f"original mutant is {original_mutant}")

    
    base_name = args.output_file[:-4] # remove .txt extension

    if args.write_pdb:
        os.makedirs(base_name + '__pdbs', exist_ok=True)
    
    logging.info(f'acceptance status = {args.peptide_acceptance}')
    logging.info(f'force status = {args.force_best_mutations}')

    annealer = TCR_pMHC_Annealer(
        hermes_dir=args.hermes_dir,
        model_version=args.model_version,
        mutant_tup=mutant_tups[0],
        peptide_L=peptide_L,
        T_regions=schedules,
        iterations=args.iters,
        region_to_optimize_energy_of=args.region_to_optimize_energy_of,
        peptide_backbone_relax=args.relax_peptide_backbone,
        peptide_acceptance=args.peptide_acceptance,
        force_best_mutations=args.force_best_mutations, 
        write_pdbdir=(base_name + '__pdbs') if args.write_pdb else None,
        use_logp_as_energy=(args.energy_to_use=='pnlogp')
    )

    from time import time
    start = time()
    annealer.anneal()
    end = time()
    print(f"Annealing took {end-start} seconds")
    
    if args.save_metric_history:
        logging.info('Saving history')
        np.save(base_name + f'.npy', annealer.metric_history)

    final_sequence = annealer.metric_history['sequences'][-2].decode('utf-8') # last one is empty for some reason

    print(annealer.metric_history['sequences'])
    print(final_sequence)

    with open(args.output_file, 'w+') as f:
        f.write(final_sequence)
    
    logging.info('Written final sampled sequence')
