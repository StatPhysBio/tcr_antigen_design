
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

    parser.add_argument('--peptide_resnum_start', type=int, default=1,
                            help='Sometimes, this is not one. This might happen for example with certain structures generated in-silico.')
    
    args = parser.parse_args()

    # read csv file
    df = pd.read_csv(args.csv_filepath)

    # load HERMES models
    trained_models_path = os.path.join(args.hermes_path, 'trained_models', args.model_version)
    model_dir_list = [os.path.join(trained_models_path, model_rel_path) for model_rel_path in os.listdir(trained_models_path)]
    models, hparams, finetuning_hparams = load_hermes_models(model_dir_list)

    # for every row in the csv file
    pnE_list = []
    pnlogp_list = []
    for i, row in df.iterrows():

        pdb = row[args.pdb_column]
        chain = row[args.chain_column]
        sequence = row[args.peptide_column]

        # use peptide mutator just to compute peptide regions
        mutator = PeptideMutator(pdb, chain, sequence, peptide_resnum_start=args.peptide_resnum_start)
        requested_regions = {'peptide': mutator.regions['peptide']} # only keep peptide region to make it quicker

        # compute pnE and pnlogp on the structure, save in new column of csv file
        ensemble_predictions_dict = predict_from_pdbfile(os.path.join(args.pdbdir, pdb + '.pdb'), models, hparams, 64, finetuning_hparams=finetuning_hparams, regions=requested_regions)

        ensemble_predictions_dict = ensemble_predictions_dict['peptide']

        ensembled_peptide_pes = np.mean(ensemble_predictions_dict['logits'], axis=0)
        if args.ensemble_at_logits_level:
            ensembled_peptide_logps = log_softmax(ensembled_peptide_pes, axis=1)
        else:
            ensembled_peptide_logps = np.log(np.mean(ensemble_predictions_dict['probabilities'], axis=0))

        ensembled_peptide_pnE = get_protein_network_energy(args.sequence, ensembled_peptide_pes, as_logp=False)
        ensembled_peptide_pnlogp = get_protein_network_energy(args.sequence, ensembled_peptide_logps, as_logp=True)

        pnE_list.append(ensembled_peptide_pnE)
        pnlogp_list.append(ensembled_peptide_pnlogp)

    # save new csv file
    df['pnE'] = np.array(pnE_list)
    df['pnlogp'] = np.array(pnlogp_list)

    df.to_csv(args.output_csv_filepath, index=False)

    
