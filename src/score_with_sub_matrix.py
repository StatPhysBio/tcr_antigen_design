
import os
import numpy as np
import pandas as pd
import argparse

from substitution_matrices import BLOSUM62

def score_with_sub_matrix(wt_peptide, mt_peptide, matrix):
    score = 0
    for wt_aa, mt_aa in zip(wt_peptide, mt_peptide):
        score += matrix[wt_aa][mt_aa]
    return score

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_csv_filepath', type=str, required=True)
    parser.add_argument('--output_csv_filepath', type=str, required=True)
    parser.add_argument('--wt_peptide_column', type=str, required=True)
    parser.add_argument('--mt_peptide_column', type=str, required=True)    
    args = parser.parse_args()

    # read csv file
    df = pd.read_csv(args.input_csv_filepath, index_col=False)

    scores = []
    for i, row in df.iterrows():
        mt_peptide = row[args.mt_peptide_column]
        wt_peptide = row[args.wt_peptide_column]

        assert len(mt_peptide) == len(wt_peptide), f"Peptides must have the same length. {mt_peptide} vs {wt_peptide} line {i}"
        
        score = score_with_sub_matrix(wt_peptide, mt_peptide, BLOSUM62)
        scores.append(score)
    
    df['substitution_matrix_score'] = np.array(scores)

    os.makedirs(os.path.dirname(args.output_csv_filepath), exist_ok=True)
    df.to_csv(args.output_csv_filepath, index=False)
