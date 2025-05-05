
import os, sys
import time
import json

import numpy as np
import pandas as pd

import logging
import argparse
from tqdm import tqdm

from pathlib import Path

from typing import *

from scipy.special import softmax, log_softmax, logsumexp

from hermes.inference.inference_hermes import predict_from_pdbfile, load_hermes_models


from get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation import PeptideMutator, get_protein_network_energy


def prediction(args, pdb, chain, sequence, models, hparams, finetuning_hparams, requested_regions=None):

    if requested_regions is None:
        # # use peptide mutator just to compute peptide regions
        # mutator = PeptideMutator(os.path.join(args.pdbdir, pdb + '.pdb'), chain, sequence, peptide_resnum_start=args.peptide_resnum_start)
        # requested_regions = {'peptide': mutator.regions['peptide']} # only keep peptide region to make it quicker
        region_ids = [(chain, resnum, ' ') for resnum in list(range(args.peptide_resnum_start, args.peptide_resnum_start + len(sequence)))]
        requested_regions = {'peptide': region_ids}

    # compute pnE and pnlogp on the structure, save in new column of csv file
    ensemble_predictions_dict = predict_from_pdbfile(os.path.join(args.pdbdir, pdb + '.pdb'), models, hparams, 64, finetuning_hparams=finetuning_hparams, regions=requested_regions)

    ensemble_predictions_dict = ensemble_predictions_dict['peptide']

    ensembled_peptide_pes = np.mean(ensemble_predictions_dict['logits'], axis=0)
    if args.ensemble_at_logits_level:
        ensembled_peptide_logps = log_softmax(ensembled_peptide_pes, axis=1)
    else:
        ensembled_peptide_logps = np.log(np.mean(ensemble_predictions_dict['probabilities'], axis=0))

    ensembled_peptide_pnE = get_protein_network_energy(sequence, ensembled_peptide_pes, as_logp=False)
    ensembled_peptide_pnlogp = get_protein_network_energy(sequence, ensembled_peptide_logps, as_logp=True)

    return ensembled_peptide_pnE, ensembled_peptide_pnlogp, requested_regions


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--hermes_path', type=str, default='/gscratch/spe/gvisan01/hermes/')
    parser.add_argument('--model_version', type=str, required=True)

    parser.add_argument('--pdbdir', type=str, required=True)
    parser.add_argument('--output_csv_filepath', type=str, required=True)

    parser.add_argument('--csv_filepath', type=str, required=True)
    parser.add_argument('--pdb_column', type=str, required=True)
    parser.add_argument('--chain_column', type=str, required=True)
    parser.add_argument('--peptide_column', type=str, required=True)

    parser.add_argument('--ensemble_at_logits_level', type=int, default=1, choices=[0, 1])

    parser.add_argument('--is_main_wt_column', type=str, default=None)

    parser.add_argument('--peptide_resnum_start', type=int, default=1,
                            help='Sometimes, this is not one. This might happen for example with certain structures generated in-silico.')
    
    parser.add_argument('--same_peptide_region', type=int, default=0, choices=[0, 1],
                        help='1 if the peptides evaluated are in the same region of the same structure; if so, toggle it to speed up the computation.')
    
    args = parser.parse_args()

    # read csv file
    if args.csv_filepath.endswith('.tsv'):
        sep = '\t'
    else:
        sep = ','
    df = pd.read_csv(args.csv_filepath, sep=sep)

    # load HERMES models
    trained_models_path = os.path.join(args.hermes_path, 'trained_models', args.model_version)
    model_dir_list = [os.path.join(trained_models_path, model_rel_path) for model_rel_path in os.listdir(trained_models_path)]
    models, hparams, finetuning_hparams = load_hermes_models(model_dir_list)

    # for every row in the csv file
    requested_regions_as_input = None
    pnE_list = []
    pnlogp_list = []
    for i, row in tqdm(df.iterrows(), total=len(df)):

        pdb = row[args.pdb_column]
        chain = row[args.chain_column].split('|')[0]
        sequence = row[args.peptide_column]

        ensembled_peptide_pnE, ensembled_peptide_pnlogp, requested_regions = prediction(args, pdb, chain, sequence, models, hparams, finetuning_hparams, requested_regions=requested_regions_as_input)

        if args.same_peptide_region:
            requested_regions_as_input = requested_regions

        pnE_list.append(ensembled_peptide_pnE)
        pnlogp_list.append(ensembled_peptide_pnlogp)

    # save new csv file
    df['pnE'] = np.array(pnE_list)
    df['pnlogp'] = np.array(pnlogp_list)

    # if args.is_main_wt_column is not None:

    #     wt_peptides = df[args.peptide_column].values[df[args.is_main_wt_column].values]

    #     pdbs_unique = np.unique(df[args.pdb_column].values)

    #     pdb_to_predictions = {}
    #     for pdb in pdbs_unique:
    #         chain = df[args.chain_column].values[df[args.pdb_column].values == pdb][0].split('|')[0]
    #         pnE_list = []
    #         pnlogp_list = []
    #         for wt_pep in wt_peptides:                    
    #             ensembled_peptide_pnE, ensembled_peptide_pnlogp = prediction(args, pdb, chain, wt_pep, models, hparams, finetuning_hparams)
    #             pnE_list.append(ensembled_peptide_pnE)
    #             pnlogp_list.append(ensembled_peptide_pnlogp)
    #         pdb_to_predictions[pdb] = (np.mean(pnE_list), np.mean(pnlogp_list))
        
    #     pne_wt = []
    #     pnlogp_wt = []
    #     for i, row in df.iterrows():
    #         pdb = row[args.pdb_column]
    #         pne_wt.append(pdb_to_predictions[pdb][0])
    #         pnlogp_wt.append(pdb_to_predictions[pdb][1])
        
    #     df[f'pnE-norm-{args.is_main_wt_column}'] = df['pnE'] - pne_wt
    #     df[f'pnlogp-norm-{args.is_main_wt_column}'] = df['pnlogp'] - pnlogp_wt

    if args.output_csv_filepath.endswith('.tsv'):
        sep = '\t'
    else:
        sep = ','
    
    df.to_csv(args.output_csv_filepath, sep=sep, index=False)

    
