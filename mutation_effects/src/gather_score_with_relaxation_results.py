
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
    parser.add_argument('--use_min_rosetta_energy_instead_of_full_average', type=int, default=0, choices=[0, 1])
    parser.add_argument('--use_min_rosetta_energy_runs_but_compute_mean', type=int, default=0, choices=[0, 1])
    args = parser.parse_args()

    assert args.use_min_rosetta_energy_instead_of_full_average + args.use_min_rosetta_energy_runs_but_compute_mean <= 1

    output_dir = os.path.join(args.experiment_dir, 'results', args.model_version)
    results_dir = os.path.join(output_dir, 'with_relaxation')
    csv_file = os.path.join(args.experiment_dir, args.csv_filename)
    output_csv_file = os.path.join(output_dir, f'{args.csv_filename.strip(".csv")}_with_relaxation-{args.model_version}-use_mt_structure=0.csv')
    
    if args.use_min_rosetta_energy_instead_of_full_average:
        results_dir = results_dir + '__runs_with_rosetta_scores'
        output_csv_file = output_csv_file.replace('_with_relaxation', '_with_relaxation_min_energy')
    elif args.use_min_rosetta_energy_runs_but_compute_mean:
        results_dir = results_dir + '__runs_with_rosetta_scores'
        output_csv_file = output_csv_file.replace('_with_relaxation', '_with_relaxation_mean_but_min_energy_runs')
    
    df = pd.read_csv(csv_file)

    # gather results, pnE and pnlogp

    for metric in ['pnE', 'pnlogp']:

        scores = []
        for i, row in df.iterrows():
            seq = row[args.seq_column]

            if args.pdbid == 'all':
                metric_files = glob(os.path.join(results_dir, f'{args.model_version}{SEP}*{SEP}{seq}{SEP}*{SEP}{metric}.npy'))
            elif args.pdbid == 'from_csv':
                pdb = row[args.pdb_column]
                metric_files = glob(os.path.join(results_dir, f'{args.model_version}{SEP}{pdb}{SEP}{seq}{SEP}*{SEP}{metric}.npy'))
            else:
                metric_files = glob(os.path.join(results_dir, f'{args.model_version}{SEP}{args.pdbid}{SEP}{seq}{SEP}*{SEP}{metric}.npy'))

            if len(metric_files) == 0:
                print(f'Warning: No metric_files found for {args.csv_filename} {args.model_version} {args.pdbid} {seq}')
                scores.append(np.nan)
                continue

            if args.use_min_rosetta_energy_instead_of_full_average:
                rosetta_energy_files = [metric_file.replace(f'{metric}', 'rosetta_energy') for metric_file in metric_files] # ensuring that the lists of files for metrics and rosetta_energies are in parallel

            temp_metric_scores = []
            temp_rosetta_energy_scores = []
            for i, metric_file in enumerate(metric_files):
                temp_metric_scores.append(np.load(metric_file))
                if args.use_min_rosetta_energy_instead_of_full_average:
                    temp_rosetta_energy_scores.append(np.load(rosetta_energy_files[i])) # assuming the lists of files for metrics and rosetta_energies are in parallel

            if args.use_min_rosetta_energy_instead_of_full_average:
                scores.append(np.hstack(temp_metric_scores)[np.argmin(np.hstack(temp_rosetta_energy_scores))])
            else:
                scores.append(np.mean(np.hstack(temp_metric_scores)))

        df[metric] = scores

    df.to_csv(output_csv_file, index=False)








    



