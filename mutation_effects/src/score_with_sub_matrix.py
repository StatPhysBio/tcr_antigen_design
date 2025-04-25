
import os
import numpy as np
import pandas as pd
import argparse

from substitution_matrices import NAME_TO_SUB_MATRIX

def score_with_sub_matrix(wt_peptide, mt_peptide, matrix):
    if 'matrix' in matrix and 'd' in matrix: # the matrix also has a position weighting!
        M = matrix['matrix']
        d = matrix['d']
        assert len(d) == len(wt_peptide), f"Length of d ({len(d)}) must match length of peptide ({len(wt_peptide)})"
    else:
        M = matrix
        d = None

    score = 0
    for i, (wt_aa, mt_aa) in enumerate(zip(wt_peptide, mt_peptide)):
        if d is not None:
            score += d[i] * M[wt_aa][mt_aa]
        else:
            score += M[wt_aa][mt_aa]

    return score

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_csv_filepath', type=str, required=True)
    parser.add_argument('--output_csv_filepath', type=str, required=True)
    parser.add_argument('--wt_peptide_column', type=str, required=True)
    parser.add_argument('--mt_peptide_column', type=str, required=True)    
    parser.add_argument('--substitution_matrix', type=str, required=True, choices=list(NAME_TO_SUB_MATRIX.keys()))  
    args = parser.parse_args()

    # read csv file
    df = pd.read_csv(args.input_csv_filepath, index_col=False)

    scores = []
    for i, row in df.iterrows():
        mt_peptide = row[args.mt_peptide_column]
        wt_peptide = row[args.wt_peptide_column]

        assert len(mt_peptide) == len(wt_peptide), f"Peptides must have the same length. {mt_peptide} vs {wt_peptide} line {i}"
        
        score = score_with_sub_matrix(wt_peptide, mt_peptide, NAME_TO_SUB_MATRIX[args.substitution_matrix])
        scores.append(score)
    
    df['substitution_matrix_score'] = np.array(scores)

    os.makedirs(os.path.dirname(args.output_csv_filepath), exist_ok=True)
    df.to_csv(args.output_csv_filepath, index=False)
