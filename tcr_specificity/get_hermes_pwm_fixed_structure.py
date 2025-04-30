
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

from hermes.utils.protein_naming import ind_to_ol_size, ol_to_ind_size


def prediction(args, pdb, chain, resnums, models, hparams, finetuning_hparams):

    region_ids = [(chain, resnum, ' ') for resnum in resnums]
    requested_regions = {'peptide': region_ids}

    # compute pnE and pnlogp on the structure, save in new column of csv file
    ensemble_predictions_dict = predict_from_pdbfile(os.path.join(args.pdbdir, pdb + '.pdb'), models, hparams, 64, finetuning_hparams=finetuning_hparams, regions=requested_regions)

    ensemble_predictions_dict = ensemble_predictions_dict['peptide']

    ensembled_peptide_pes = np.mean(ensemble_predictions_dict['logits'], axis=0)
    if args.ensemble_at_logits_level:
        ensembled_peptide_probas = softmax(ensembled_peptide_pes, axis=1)
    else:
        ensembled_peptide_probas = np.mean(ensemble_predictions_dict['probabilities'], axis=0)

    return ensembled_peptide_probas


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--hermes_path', type=str, default='/gscratch/spe/gvisan01/hermes/')
    parser.add_argument('--model_version', type=str, default='hermes_py_000')

    parser.add_argument('--pdbdir', type=str, default='pdbs/pmhc')
    parser.add_argument('--output_dir', type=str, default='pwm_csv_files/mhc_hcnn_fixed_structure')

    parser.add_argument('--input_csv_filepath', type=str, default='pdbs_allele_df_with_peptide_info.csv')
    parser.add_argument('--pdb_column', type=str, default='pdbid')
    parser.add_argument('--chain_column', type=str, default='peptide_chain')
    parser.add_argument('--resnums_column', type=str, default='peptide_resnums')
    parser.add_argument('--sequence_column', type=str, default='wt_peptide')
    parser.add_argument('--allele_column', type=str, default='mhc_allele')

    parser.add_argument('--ensemble_at_logits_level', type=int, default=1, choices=[0, 1])

    parser.add_argument('--peptide_resnum_start', type=int, default=1,
                            help='Sometimes, this is not one. This might happen for example with certain structures generated in-silico.')

    args = parser.parse_args()

    output_dir = os.path.join(args.output_dir, args.model_version)
    os.makedirs(output_dir, exist_ok=True)

    # read csv file
    df = pd.read_csv(args.input_csv_filepath)

    # load HERMES models
    trained_models_path = os.path.join(args.hermes_path, 'trained_models', args.model_version)
    model_dir_list = [os.path.join(trained_models_path, model_rel_path) for model_rel_path in os.listdir(trained_models_path)]
    models, hparams, finetuning_hparams = load_hermes_models(model_dir_list)

    # for every row in the csv file
    for i, row in tqdm(df.iterrows(), total=len(df)):

        pdb = row[args.pdb_column]
        chain = row[args.chain_column]
        resnums = map(int, row[args.resnums_column].split('|'))
        wt_peptide = row[args.sequence_column]
        allele = row[args.allele_column] 

        peptide_pwm = prediction(args, pdb, chain, resnums, models, hparams, finetuning_hparams)

        filename = f'{pdb}__{wt_peptide}__{allele}.csv'
        filepath = os.path.join(output_dir, filename)

        df = pd.DataFrame(peptide_pwm, index=range(1, peptide_pwm.shape[0]+1), columns=[ind_to_ol_size[ind] for ind in range(20)])

        df.to_csv(filepath)


