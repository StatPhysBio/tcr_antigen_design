
import os
from glob import glob
import numpy as np
import pandas as pd
from tqdm import tqdm

import argparse

# from hermes.utils.protein_naming import ol_to_ind_size, ind_to_ol_size
ind_to_ol_size = {0: 'G', 1: 'A', 2: 'C', 3: 'S', 4: 'P', 5: 'T', 6: 'V', 7: 'D', 8: 'I', 9: 'L', 10: 'N', 11: 'M', 12: 'Q', 13: 'K', 14: 'E', 15: 'H', 16: 'F', 17: 'R', 18: 'Y', 19: 'W'}
ol_to_ind_size = {v: k for k, v in ind_to_ol_size.items()}

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

SEP = '$'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--csv_filename', type=str, required=True)
    parser.add_argument('--pdbid', type=str, required=True)
    parser.add_argument('--model_version', type=str, required=True, help='hermes model version')
    parser.add_argument('--pdb_column', type=str, default='pdb', help='Only if pdbid is from_csv')
    parser.add_argument('--seq_column', type=str, default='sequence', help='column name of the sequence in the csv file')
    parser.add_argument('--output_in_inputfile', type=int, default=0, choices=[0, 1])
    parser.add_argument('--use_max_instead_of_mean', type=int, default=0, choices=[0, 1])
    parser.add_argument('--use_min_rosetta_energy_instead_of_full_average', type=int, default=0, choices=[0, 1])
    parser.add_argument('--use_min_rosetta_energy_runs_but_compute_mean', type=int, default=0, choices=[0, 1])
    args = parser.parse_args()

    assert args.use_min_rosetta_energy_instead_of_full_average + args.use_min_rosetta_energy_runs_but_compute_mean <= 1

    output_dir = os.path.join('hermes_scores', args.model_version)
    results_dir = os.path.join(output_dir, 'with_relaxation')

    file_extension = args.csv_filename[-4:]
    if file_extension == '.csv':
        sep = ','
    elif file_extension == '.tsv':
        sep = '\t'
    
    if args.output_in_inputfile:
        output_csv_file = args.csv_filename
        output_dir = os.path.join(os.path.dirname(output_csv_file), args.model_version)
        results_dir = os.path.join(output_dir, 'with_relaxation')
    else:
        output_csv_file = os.path.join(output_dir, f'{args.csv_filename[:-4]}__{args.model_version}__relaxed{file_extension}')
    
    if args.use_min_rosetta_energy_instead_of_full_average:
        results_dir = results_dir + '__runs_with_rosetta_scores'
        output_csv_file = output_csv_file.replace('_with_relaxation', '_with_relaxation_min_energy')
    elif args.use_min_rosetta_energy_runs_but_compute_mean:
        results_dir = results_dir + '__runs_with_rosetta_scores'
        output_csv_file = output_csv_file.replace('_with_relaxation', '_with_relaxation_mean_but_min_energy_runs')
    
    df = pd.read_csv(args.csv_filename, sep=sep)

    # gather results, pnE and pnlogp

    for metric in ['pnE', 'pnlogp']:

        scores = []
        for i, row in tqdm(df.iterrows(), total=len(df)):
            seq = row[args.seq_column]

            if args.pdbid == 'all':
                metric_files = glob(os.path.join(results_dir, f'{args.model_version}{SEP}*{SEP}{seq}{SEP}*{SEP}{metric}.npy'))
            elif args.pdbid == 'from_csv':
                pdb = row[args.pdb_column]
                metric_files = glob(os.path.join(results_dir, f'{args.model_version}{SEP}{pdb}{SEP}{seq}{SEP}*{SEP}{metric}.npy'))
            else:
                metric_files = glob(os.path.join(results_dir, f'{args.model_version}{SEP}{args.pdbid}{SEP}{seq}{SEP}*{SEP}{metric}.npy'))

            if len(metric_files) == 0:
                print(f'Warning: No metric_files found for {args.model_version} {args.pdbid} {seq}')
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
                if args.use_max_instead_of_mean:
                    scores.append(np.max(np.hstack(temp_metric_scores)))
                else:
                    scores.append(np.mean(np.hstack(temp_metric_scores)))

        df[metric] = scores

    df.to_csv(output_csv_file, sep=sep, index=False)








    



