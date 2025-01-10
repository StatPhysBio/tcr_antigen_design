

import os
import pandas as pd
from tqdm import tqdm

from hermes.utils.protein_naming import ind_to_ol_size, ol_to_ind_size

from constants import WT_SEQ

# from Bio.Align import substitution_matrices
# BLOSUM62 = substitution_matrices.load('BLOSUM62')

BLOSUM62_DICT = {
    'C': {'C': 9, 'S': -1, 'T': -1, 'A': 0, 'G': -3, 'P': -3, 'D': -3, 'E': -4, 'Q': -3, 'N': -3, 'H': -3, 'R': -3, 'K': -3, 'M': -1, 'I': -1, 'L': -1, 'V': -1, 'W': -2, 'Y': -2, 'F': -2},
    'S': {'C': -1, 'S': 4, 'T': 1, 'A': 1, 'G': 0, 'P': -1, 'D': 0, 'E': 0, 'Q': 0, 'N': 1, 'H': -1, 'R': -1, 'K': 0, 'M': -1, 'I': -2, 'L': -2, 'V': -2, 'W': -3, 'Y': -2, 'F': -2},
    'T': {'C': -1, 'S': 1, 'T': 5, 'A': 0, 'G': -2, 'P': -1, 'D': -1, 'E': -1, 'Q': -1, 'N': 0, 'H': -2, 'R': -1, 'K': -1, 'M': -1, 'I': -1, 'L': -1, 'V': 0, 'W': -2, 'Y': -2, 'F': -2},
    'A': {'C': 0, 'S': 1, 'T': 0, 'A': 4, 'G': 0, 'P': -1, 'D': -2, 'E': -1, 'Q': -1, 'N': -2, 'H': -2, 'R': -1, 'K': -1, 'M': -1, 'I': -1, 'L': -1, 'V': 0, 'W': -3, 'Y': -2, 'F': -2},
    'G': {'C': -3, 'S': 0, 'T': -2, 'A': 0, 'G': 6, 'P': -2, 'D': -1, 'E': -2, 'Q': -2, 'N': 0, 'H': -2, 'R': -2, 'K': -2, 'M': -3, 'I': -4, 'L': -4, 'V': -3, 'W': -2, 'Y': -3, 'F': -3},
    'P': {'C': -3, 'S': -1, 'T': -1, 'A': -1, 'G': -2, 'P': 7, 'D': -1, 'E': -1, 'Q': -1, 'N': -1, 'H': -2, 'R': -2, 'K': -1, 'M': -2, 'I': -3, 'L': -3, 'V': -2, 'W': -4, 'Y': -3, 'F': -4},
    'D': {'C': -3, 'S': 0, 'T': -1, 'A': -2, 'G': -1, 'P': -1, 'D': 6, 'E': 2, 'Q': 0, 'N': 1, 'H': -1, 'R': -2, 'K': -1, 'M': -3, 'I': -3, 'L': -4, 'V': -3, 'W': -4, 'Y': -3, 'F': -3},
    'E': {'C': -4, 'S': 0, 'T': -1, 'A': -1, 'G': -2, 'P': -1, 'D': 2, 'E': 5, 'Q': 2, 'N': 0, 'H': 0, 'R': 0, 'K': 1, 'M': -2, 'I': -3, 'L': -3, 'V': -2, 'W': -3, 'Y': -2, 'F': -3},
    'Q': {'C': -3, 'S': 0, 'T': -1, 'A': -1, 'G': -2, 'P': -1, 'D': 0, 'E': 2, 'Q': 5, 'N': 0, 'H': 0, 'R': 1, 'K': 1, 'M': 0, 'I': -3, 'L': -2, 'V': -2, 'W': -2, 'Y': -1, 'F': -3},
    'N': {'C': -3, 'S': 1, 'T': 0, 'A': -2, 'G': 0, 'P': -2, 'D': 1, 'E': 0, 'Q': 0, 'N': 6, 'H': 1, 'R': 0, 'K': 0, 'M': -2, 'I': -3, 'L': -3, 'V': -3, 'W': -4, 'Y': -2, 'F': -3},
    'H': {'C': -3, 'S': -1, 'T': -2, 'A': -2, 'G': -2, 'P': -2, 'D': -1, 'E': 0, 'Q': 0, 'N': 1, 'H': 8, 'R': 0, 'K': -1, 'M': -2, 'I': -3, 'L': -3, 'V': -3, 'W': -2, 'Y': 2, 'F': -1},
    'R': {'C': -3, 'S': -1, 'T': -1, 'A': -1, 'G': -2, 'P': -2, 'D': -2, 'E': 0, 'Q': 1, 'N': 0, 'H': 0, 'R': 5, 'K': 2, 'M': -1, 'I': -3, 'L': -2, 'V': -3, 'W': -3, 'Y': -2, 'F': -3},
    'K': {'C': -3, 'S': 0, 'T': -1, 'A': -1, 'G': -2, 'P': -1, 'D': -1, 'E': 1, 'Q': 1, 'N': 0, 'H': -1, 'R': 2, 'K': 5, 'M': -1, 'I': -3, 'L': -2, 'V': -2, 'W': -3, 'Y': -2, 'F': -3},
    'M': {'C': -1, 'S': -1, 'T': -1, 'A': -1, 'G': -3, 'P': -2, 'D': -3, 'E': -2, 'Q': 0, 'N': -2, 'H': -2, 'R': -1, 'K': -1, 'M': 5, 'I': 1, 'L': 2, 'V': 1, 'W': -1, 'Y': -1, 'F': 0},
    'I': {'C': -1, 'S': -2, 'T': -1, 'A': -1, 'G': -4, 'P': -3, 'D': -3, 'E': -3, 'Q': -3, 'N': -3, 'H': -3, 'R': -3, 'K': -3, 'M': 1, 'I': 4, 'L': 2, 'V': 3, 'W': -3, 'Y': -1, 'F': 0},
    'L': {'C': -1, 'S': -2, 'T': -1, 'A': -1, 'G': -4, 'P': -3, 'D': -4, 'E': -3, 'Q': -2, 'N': -3, 'H': -3, 'R': -2, 'K': -2, 'M': 2, 'I': 2, 'L': 4, 'V': 1, 'W': -2, 'Y': -1, 'F': 0},
    'V': {'C': -1, 'S': -2, 'T': 0, 'A': 0, 'G': -3, 'P': -2, 'D': -3, 'E': -2, 'Q': -2, 'N': -3, 'H': -3, 'R': -3, 'K': -2, 'M': 1, 'I': 3, 'L': 1, 'V': 4, 'W': -3, 'Y': -1, 'F': -1},
    'W': {'C': -2, 'S': -3, 'T': -2, 'A': -3, 'G': -2, 'P': -4, 'D': -4, 'E': -3, 'Q': -2, 'N': -4, 'H': -2, 'R': -3, 'K': -3, 'M': -1, 'I': -3, 'L': -2, 'V': -3, 'W': 11, 'Y': 2, 'F': 1},
    'Y': {'C': -2, 'S': -2, 'T': -2, 'A': -2, 'G': -3, 'P': -3, 'D': -3, 'E': -2, 'Q': -1, 'N': -2, 'H': 2, 'R': -2, 'K': -2, 'M': -1, 'I': -1, 'L': -1, 'V': -1, 'W': 2, 'Y': 7, 'F': 3},
    'F': {'C': -2, 'S': -2, 'T': -2, 'A': -2, 'G': -3, 'P': -4, 'D': -3, 'E': -3, 'Q': -3, 'N': -3, 'H': -1, 'R': -3, 'K': -3, 'M': 0, 'I': 0, 'L': 0, 'V': -1, 'W': 1, 'Y': 3, 'F': 6}
}

def augment_df_with_blosum(df, wt_seq, seq_column='peptide'):
    df['blosum_full'] = [sum([BLOSUM62_DICT[wt_seq[i]][seq[i]] for i in range(len(seq))]) for seq in df[seq_column]]
    df['blosum_diff'] = [sum([BLOSUM62_DICT[wt_seq[i]][seq[i]] for i in range(len(seq)) if seq[i] != wt_seq[i]]) / sum([1 for i in range(len(seq)) if seq[i] != wt_seq[i]]) if seq != wt_seq else 0.0 for seq in df[seq_column]]
    return df

def combine_scores(df, col1, col2, rescale='min_max', weight_of_col1=0.5):
    combined_score_name = f'{col1}__{col2}__{rescale}__{weight_of_col1}'
    if rescale == 'min_max':
        df[combined_score_name] = (df[col1] - df[col1].min()) / ((df[col1].max() - df[col1].min())) * weight_of_col1 + ((df[col2] - df[col2].min()) / (df[col2].max() - df[col2].min())) * (1 - weight_of_col1)
    elif rescale == 'z_score':
        df[combined_score_name] = ((df[col1] - df[col1].mean()) / df[col1].std()) * weight_of_col1 + ((df[col2] - df[col2].mean()) / df[col2].std()) * (1 - weight_of_col1)
    return df


if __name__ == '__main__':

    tsvfiles = ['wildtype/wildtype_w_pae.tsv',
                'blosum62/sample_peptides_from_blosum62__temperature=1.0_w_pae.tsv',
                'blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae.tsv',
                'blosum62/sample_peptides_from_blosum62__temperature=3.0_w_pae.tsv',
                'hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae.tsv',
                'hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae.tsv',
                'hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae.tsv',
                'hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae.tsv',
                'mhc_pwm/mhc_motif_peptides_w_pae.tsv',
                'proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7_w_pae.tsv',
                'proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7_w_pae.tsv']

    for tsvfile in tqdm(tsvfiles):
        df = pd.read_csv(tsvfile, sep='\t')

        ## make neg_pmhc_tcr_pae
        df['neg_pmhc_tcr_pae'] = -df['pmhc_tcr_pae']

        ## augment dataframes with blosum scores
        df = augment_df_with_blosum(df, WT_SEQ)

        ## combine scores
        df = combine_scores(df, 'neg_pmhc_tcr_pae', 'blosum_full', rescale='min_max', weight_of_col1=0.5)
        df = combine_scores(df, 'neg_pmhc_tcr_pae', 'blosum_diff', rescale='min_max', weight_of_col1=0.5)

        ## save the augmented dataframe
        df.to_csv(tsvfile.replace('.tsv', '_w_blosum.tsv'), sep='\t', index=False)


