
import os
from glob import glob
import numpy as np
import pandas as pd

import argparse

from hermes.utils.protein_naming import ol_to_ind_size, ind_to_ol_size

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

SEP = '$'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--experiment_dir', type=str, required=True)
    parser.add_argument('--csv_filename', type=str, required=True)
    parser.add_argument('--pdbid', type=str, required=True)
    parser.add_argument('--model_version', type=str, required=True, help='hermes model version')
    parser.add_argument('--pdb_column', type=str, default='pdb', help='Only if pdbid is from_csv')
    parser.add_argument('--seq_column', type=str, default='sequence', help='column name of the sequence in the csv file')
    args = parser.parse_args()

    output_dir = os.path.join(args.experiment_dir, 'results', args.model_version)
    results_dir = os.path.join(output_dir, 'with_relaxation')
    csv_file = os.path.join(args.experiment_dir, args.csv_filename)
    output_csv_file = os.path.join(output_dir, f'{args.csv_filename.strip(".csv")}_with_relaxation-{args.model_version}-use_mt_structure=0.csv')

    df = pd.read_csv(csv_file)

    # gather results, pnE ad pnlogp

    for metric in ['pnE', 'pnlogp']:

        scores = []
        for i, row in df.iterrows():
            seq = row[args.seq_column]

            if args.pdbid == 'all':
                files = glob(os.path.join(results_dir, f'{args.model_version}{SEP}*{SEP}{seq}{SEP}*{SEP}{metric}.npy'))
            elif args.pdbid == 'from_csv':
                pdb = row[args.pdb_column]
                files = glob(os.path.join(results_dir, f'{args.model_version}{SEP}{pdb}{SEP}{seq}{SEP}*{SEP}{metric}.npy'))
            else:
                files = glob(os.path.join(results_dir, f'{args.model_version}{SEP}{args.pdbid}{SEP}{seq}{SEP}*{SEP}{metric}.npy'))

            assert len(files) > 0, f'No files found for {args.model_version} {args.pdbid} {seq}'

            temp_scores = []
            for file in files:
                temp_scores.append(np.load(file))
            scores.append(np.mean(np.hstack(temp_scores)))

        df[metric] = scores

    df.to_csv(output_csv_file, index=False)








    



