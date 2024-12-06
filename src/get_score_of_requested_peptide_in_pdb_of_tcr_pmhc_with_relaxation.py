

import os, sys
import time
import json

import numpy as np

import logging
import argparse
from tqdm import tqdm

import pyrosetta

from pathlib import Path

from typing import *

from scipy.special import softmax, log_softmax, logsumexp

from hermes.inference.inference_hermes import predict_from_pdbfile, load_hermes_models

from hermes.utils.protein_naming import ind_to_ol_size, ol_to_ind_size

import pyrosetta_utils as pru

from pyrosetta_utils import (
    default_flags,
    fastrelax_positions, make_mutations)

pyrosetta.init(default_flags, silent=True)

from pocket import (
    get_tcr_pmhc_peptide_pocket_resnums, get_pockets)

from pep import pep_to_pdb



class PeptideMutator(object):
    def __init__(self,
                 pdbpath,
                 chain,
                 sequence,
                 relax_peptide_backbone=True,
                 peptide_resnum_start=1
                 ):
        super().__init__()
    
        self.init_score_fxns()
    
        # prepare protein pose    
        self.pdbpath = os.path.join(args.pdbdir, args.pdb+'.pdb')
        self.pose = pyrosetta.pose_from_pdb(self.pdbpath)


        # get peptide residue ids
        peptide_seqnums = []
        for i in range(peptide_resnum_start, peptide_resnum_start + len(sequence)):
            seqnum = pru.get_pose_residue_number(self.pose, chain, i)
            if seqnum != 0:
                peptide_seqnums.append(seqnum)
        
        ids = [pru.get_pdb_residue_info(self.pose,i) for i in peptide_seqnums]
        print(ids)
        
        N_mutants = 1
        peptide_L = len(sequence)
        
        if sequence == "":
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
        self.mutant_tup = mutant_tups[0]

        self.relax_peptide_backbone = relax_peptide_backbone

        # extract info from mutant tup
        self.pdbfile = self.mutant_tup[0]
        self.pdb = Path(self.pdbfile).stem
        self.chain = self.mutant_tup[1]
        self.mutation_inds = self.mutant_tup[2]    
        self.init_seq = ''.join(self.mutant_tup[3])
        self.mutant_tup = self.mutant_tup

        print('chain:', self.chain)

        self.init_original_pose()
        
    
    def init_original_pose(self):
        self.pose = pyrosetta.pose_from_pdb(self.pdbpath)

        self.peptide_resnums, self.pocket_resnums, self.complex_resnums = (
            get_tcr_pmhc_peptide_pocket_resnums(
                self.pose, self.mutation_inds, self.chain))
        self.pocket_size = len(self.pocket_resnums)

        self.peptide_resnums = list(self.peptide_resnums)
        self.complex_resnums = list(self.complex_resnums)
        self.pocket_resnums = list(self.pocket_resnums)

        self.regions = get_pockets(self.pdb, None, self.peptide_resnums, self.pose) # None -> self.pocket_resnums
    
    def init_score_fxns(self):
        # scorefunctions
        # self.scorefxn = pyrosetta.create_score_function('ref2015.wts') # --> never using this
        self.scorefxn_cartesian = pyrosetta.create_score_function('ref2015_cart.wts')
    
    def mutate(self, verbose=False):
        mutations = {
            seqpos : mutation for seqpos,mutation in zip(self.peptide_resnums, self.mutant_tup[3])
        }
        logging.info(f"Making mutations to {self.mutant_tup[3]}")
        make_mutations(mutations, self.pose, verbose=verbose)
        logging.info("relaxing structure")
        logging.debug(f"Backbone relax sites are {self.peptide_resnums * self.relax_peptide_backbone}")
        print(self.peptide_resnums * self.relax_peptide_backbone)
        print(self.complex_resnums)
        fastrelax_positions(
            self.scorefxn_cartesian,
            self.peptide_resnums * self.relax_peptide_backbone,
            self.complex_resnums, # <-- complex_resnums is the full pocket within 10. AA of the peptide
            self.pose,
            nrepeats = 1
        )


def collect_pnEs_for_sequences(pne_dir, sequences):
    sequences = set(sequences)
    all_files = os.listdir(pne_dir)
    pnE_dict = {}
    for file in all_files:
        if 'pnE_ensemble' not in file:
            continue
        split_file = file.split('_')
        pdb, sequence = split_file[0], split_file[1]
        if sequence in sequences:
            if sequence not in pnE_dict:
                pnE_dict[sequence] = {}
            pnE_dict[sequence][pdb] = np.load(os.path.join(pne_dir, file))
    return pnE_dict

def collect_pess_for_sequences(pes_dir, sequences):
    sequences = set(sequences)
    all_files = os.listdir(pes_dir)
    pes_dict = {}
    for file in all_files:
        if 'pes_ensemble' not in file:
            continue
        split_file = file.split('_')
        pdb, sequence = split_file[0], split_file[1]
        if sequence in sequences:
            if sequence not in pes_dict:
                pes_dict[sequence] = {}
            pes_dict[sequence][pdb] = np.load(os.path.join(pes_dir, file))
    return pes_dict

def collect_logpss_for_sequences(pes_dir, sequences):
    sequences = set(sequences)
    all_files = os.listdir(pes_dir)
    pes_dict = {}
    for file in all_files:
        if 'logps_ensemble' not in file:
            continue
        split_file = file.split('_')
        pdb, sequence = split_file[0], split_file[1]
        if sequence in sequences:
            if sequence not in pes_dict:
                pes_dict[sequence] = {}
            pes_dict[sequence][pdb] = np.load(os.path.join(pes_dir, file))
    return pes_dict

def get_mean_and_std_pnE_per_sequence(pnE_dict):
    temp_pnE_dict = {}
    mean_pnE_dict = {}
    std_pnE_dict = {}
    for sequence in pnE_dict:
        temp_pnE_dict[sequence] = []
        for pdb in pnE_dict[sequence]:
            temp_pnE_dict[sequence].append(pnE_dict[sequence][pdb])
        mean_pnE_dict[sequence] = np.mean(np.hstack(temp_pnE_dict[sequence]))
        std_pnE_dict[sequence] = np.std(np.hstack(temp_pnE_dict[sequence]))
    return mean_pnE_dict, std_pnE_dict

def get_mean_and_std_pnE_per_sequence_for_pdbs(pnE_dict, allowed_pdbs):
    temp_pnE_dict = {}
    mean_pnE_dict = {}
    std_pnE_dict = {}
    for sequence in pnE_dict:
        temp_pnE_dict[sequence] = []
        for pdb in pnE_dict[sequence]:
            if pdb in allowed_pdbs:
                temp_pnE_dict[sequence].append(pnE_dict[sequence][pdb])
        mean_pnE_dict[sequence] = np.mean(np.hstack(temp_pnE_dict[sequence]))
        std_pnE_dict[sequence] = np.std(np.hstack(temp_pnE_dict[sequence]))
    return mean_pnE_dict, std_pnE_dict
            
def get_mean_and_std_pnE_per_sequence_for_sequence_pdb_pairs(pnE_dict, allowed_sequence_pdb_dict):
    temp_pnE_dict = {}
    mean_pnE_dict = {}
    std_pnE_dict = {}
    for sequence in pnE_dict:
        if sequence in allowed_sequence_pdb_dict:
            temp_pnE_dict[sequence] = []
            for pdb in pnE_dict[sequence]:
                if pdb in allowed_sequence_pdb_dict[sequence]:
                    temp_pnE_dict[sequence].append(pnE_dict[sequence][pdb])
            mean_pnE_dict[sequence] = np.mean(np.hstack(temp_pnE_dict[sequence]))
            std_pnE_dict[sequence] = np.std(np.hstack(temp_pnE_dict[sequence]))
    return mean_pnE_dict, std_pnE_dict


def get_mean_logpmt_minus_logpwt_per_sequence(pes_dict, wt_sequence, get_pseudoenergy_instead=False):

    def get_indices_of_variance(pmt, pwt):
        indices = []
        for i, (mt, wt) in enumerate(zip(pmt, pwt)):
            if mt != wt:
                indices.append(i)
        if len(indices) == 0: return None
        else: return np.array(indices)
    
    # get logp of wt
    # wt_pes = np.mean(np.concatenate([pes_dict[wt_sequence][pdb] for pdb in pes_dict[wt_sequence]], axis=0), axis=0)
    wt_pes = np.mean(pes_dict[wt_sequence][pep_to_pdb[wt_sequence]], axis=0)
    if get_pseudoenergy_instead:
        wt_logp = wt_pes
    else:
        wt_logp = log_softmax(wt_pes, axis=1)
    wt_indices = np.array([ol_to_ind_size[ol] for ol in wt_sequence])
    wt_logp = wt_logp[np.arange(len(wt_indices)), wt_indices]
    
    # get logp of mutants
    temp_dict = {}
    logpmt_minus_logpwt_dict = {}
    for sequence in pes_dict:
        temp_dict[sequence] = []
        # for pdb in pes_dict[sequence]:
        #     temp_dict[sequence].append(pes_dict[sequence][pdb])
        temp_dict[sequence].append(pes_dict[sequence][pep_to_pdb[wt_sequence]])
        temp_dict[sequence] = np.mean(np.concatenate(temp_dict[sequence], axis=0), axis=0)

        indices_of_variance = get_indices_of_variance(sequence, wt_sequence)

        if indices_of_variance is None:
            logpmt_minus_logpwt_dict[sequence] = 0.0 # wildtype!
        # elif len(indices_of_variance) > 1 or len(indices_of_variance) == 0:
        #     logpmt_minus_logpwt_dict[sequence] = np.nan
        else:
            indices = np.array([ol_to_ind_size[ol] for ol in sequence])
            if get_pseudoenergy_instead:
                logpmt = temp_dict[sequence]
            else:
                logpmt = np.log(softmax(temp_dict[sequence], axis=1))
            logpmt = logpmt[np.arange(len(indices)), indices]
            logpmt_minus_logpwt_dict[sequence] = np.mean(logpmt[indices_of_variance] - wt_logp[indices_of_variance])
    
    return logpmt_minus_logpwt_dict


def get_protein_network_energy(
    seq: str,
    pseudoenergies: np.ndarray,
    as_logp: bool = False
) -> float:
    """
    Get protein network energy for a given sequence over a matrix of 
    predictions
    """
    idx_seq = [ol_to_ind_size[x] for x in seq]
    if as_logp:
        return np.sum((pseudoenergies - logsumexp(pseudoenergies, axis=1, keepdims=True))[np.arange(len(seq)),idx_seq])
    else:
        return np.sum(pseudoenergies[np.arange(len(seq)),idx_seq])


if __name__ == '__main__':

    # get pdb and sequence I wanna put into the peptide (which is always at chain C

    parser = argparse.ArgumentParser()
    parser.add_argument('--hermes_path', type=str, default='/gscratch/spe/gvisan01/hermes/')
    parser.add_argument('--model_version', type=str, required=True)
    
    parser.add_argument('--pdb', type=str, required=True)
    parser.add_argument('--pdbdir', type=str, default='/gscratch/spe/gvisan01/tcr_pmhc/pdbs')

    parser.add_argument('--chain', type=str, required=True,
                        help='The chain corresponding to the peptide within the pdb. Most commonly C.')
    
    parser.add_argument('--sequence', type=str, required=True,
                        help='Desired peptide sequence to mutate into the PDB')
    
    parser.add_argument('--peptide_resnum_start', type=int, default=1,
                            help='Sometimes, this is not one. This might happen for example with certain structures generated in-silico.')

    parser.add_argument('--ensemble_at_logits_level', type=int, default=1, choices=[0, 1])

    parser.add_argument('--output_dir', type=str, default='./peptide_pnEs')

    parser.add_argument('--relax_peptide_backbone', type=int, default=1, choices=[0, 1]) 

    parser.add_argument('--num_repeats', type=int, default=1)

    parser.add_argument('--job', type=str, default=0)

    parser.add_argument('--divisor_for_name', type=str, default='$')

    parser.add_argument('--verbose', type=int, default=0, choices=[0, 1])

    parser.add_argument('--logging', dest='logging', default='DEBUG', type=str)
    parser.add_argument('--pyrosetta_logging', dest='pyrosetta_logging', default=None)
    parser.add_argument('--rosetta_logging', dest='rosetta_logging', default=None)
    args = parser.parse_args()

    logging.getLogger().setLevel(args.logging)
    if args.pyrosetta_logging != None:
        logging.getLogger("pyrosetta").setLevel(args.pyrosetta_logging)
    if args.rosetta_logging != None:
        logging.getLogger("rosetta").setLevel(args.rosetta_logging)
    
    # make mutator object
    mutator = PeptideMutator(args.pdb, args.chain, args.sequence, relax_peptide_backbone=args.relax_peptide_backbone, peptide_resnum_start=args.peptide_resnum_start)

    # load HERMES models
    trained_models_path = os.path.join(args.hermes_path, 'trained_models', args.model_version)
    model_dir_list = [os.path.join(trained_models_path, model_rel_path) for model_rel_path in os.listdir(trained_models_path)]
    models, hparams, finetuning_hparams = load_hermes_models(model_dir_list)

    pnE_ensemble = []
    pnlogp_ensemble = []
    pes_ensemble = []
    logps_ensemble = []
    for i in tqdm(range(args.num_repeats)):
        mutator.init_original_pose()
        mutator.mutate(args.verbose)

        requested_regions = {'peptide': mutator.regions['peptide']} # only keep peptide region to make it quicker
        print(requested_regions['peptide'])

        # wrie temp pdbfile with current relaxation (easiest way with current code, not fastest way)
        import tempfile
        with tempfile.NamedTemporaryFile(delete=True, suffix='.pdb') as f:
            mutator.pose.dump_pdb(f.name)
            ensemble_predictions_dict = predict_from_pdbfile(f.name, models, hparams, 64, finetuning_hparams=finetuning_hparams, regions=requested_regions)
        
        # tempfile_name = 'tempfile.pdb'
        # mutator.pose.dump_pdb(tempfile_name)
        # ensemble_predictions_dict = predict_from_pdbfile(tempfile_name, models, hparams, 64, finetuning_hparams=finetuning_hparams, regions=requested_regions)
        
        ensemble_predictions_dict = ensemble_predictions_dict['peptide']

        # print(ensemble_predictions_dict['logits'])
        # print(ensemble_predictions_dict['logits'].shape)

        # print(args.sequence)
        # print(len(args.sequence))

        ensembled_peptide_pes = np.mean(ensemble_predictions_dict['logits'], axis=0)
        if args.ensemble_at_logits_level:
            ensembled_peptide_logps = log_softmax(ensembled_peptide_pes, axis=1)
        else:
            ensembled_peptide_logps = np.log(np.mean(ensemble_predictions_dict['probabilities'], axis=0))

        ensembled_peptide_pnE = get_protein_network_energy(args.sequence, ensembled_peptide_pes, as_logp=False)
        ensembled_peptide_pnlogp = get_protein_network_energy(args.sequence, ensembled_peptide_logps, as_logp=True)

        pnE_ensemble.append(ensembled_peptide_pnE)
        pnlogp_ensemble.append(ensembled_peptide_pnlogp)
        pes_ensemble.append(ensembled_peptide_pes)
        logps_ensemble.append(ensembled_peptide_logps)

    pnE_ensemble = np.array(pnE_ensemble)
    pnlogp_ensemble = np.array(pnlogp_ensemble)
    pes_ensemble = np.stack(pes_ensemble, axis=0)
    logps_ensemble = np.stack(logps_ensemble, axis=0)

    os.makedirs(args.output_dir, exist_ok=True)

    DIV = args.divisor_for_name
    
    # save results
    os.makedirs(args.output_dir, exist_ok=True)
    np.save(os.path.join(args.output_dir, f'{args.model_version}{DIV}{args.pdb}{DIV}{args.sequence}{DIV}{args.job}{DIV}{args.num_repeats}{DIV}pnE.npy'), pnE_ensemble)
    np.save(os.path.join(args.output_dir, f'{args.model_version}{DIV}{args.pdb}{DIV}{args.sequence}{DIV}{args.job}{DIV}{args.num_repeats}{DIV}pnlogp.npy'), pnlogp_ensemble)
    # np.save(os.path.join(args.output_dir, f'{args.model_version}{DIV}{args.pdb}{DIV}{args.sequence}{DIV}{args.job}{DIV}{args.num_repeats}{DIV}pes.npy'), pes_ensemble)
    # np.save(os.path.join(args.output_dir, f'{args.model_version}{DIV}{args.pdb}{DIV}{args.sequence}{DIV}{args.job}{DIV}{args.num_repeats}{DIV}logps.npy'), logps_ensemble)

