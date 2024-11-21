

import os
import pandas as pd
from tqdm import tqdm

from protein_holography_pytorch.utils.protein_naming import ind_to_ol_size, ol_to_ind_size

from constants import WT_SEQ

from Bio.Align import substitution_matrices

BLOSUM62 = substitution_matrices.load('BLOSUM62')

def augment_df_with_blosum(df, wt_seq, seq_column='peptide'):
    df['blosum_full'] = [sum([BLOSUM62[wt_seq[i], seq[i]] for i in range(len(seq))]) for seq in df[seq_column]]
    df['blosum_diff'] = [sum([BLOSUM62[wt_seq[i], seq[i]] for i in range(len(seq)) if seq[i] != wt_seq[i]]) / sum([1 for i in range(len(seq)) if seq[i] != wt_seq[i]]) if seq != wt_seq else 0.0 for seq in df[seq_column]]
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
                # 'blosum62/sample_peptides_from_blosum62__temperature=1.0_w_pae.tsv',
                # 'blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae.tsv',
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


