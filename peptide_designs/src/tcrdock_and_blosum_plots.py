

import os, sys
import numpy as np
import pandas as pd

from scipy.stats import pearsonr, spearmanr

import matplotlib.pyplot as plt

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# from protein_holography_pytorch.utils.protein_naming import ind_to_ol_size, ol_to_ind_size
from hermes.utils.protein_naming import ind_to_ol_size, ol_to_ind_size

import argparse


'''

From tcrdock colab notebook: "Models with PAE values less than 6.5 or 7 are higher confidence; models with PAE values greater than 7.5 or 8 are low confidence."

'''


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--wt_seq', type=str, required=True, help='wildtype peptide sequence')
    args = parser.parse_args()

    WT_SEQ = args.wt_seq


    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan', 'blue', 'orange', 'green', 'red', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']

    df_wildtype = pd.read_csv('wildtype/wildtype_w_pae_w_blosum.tsv', sep='\t')

    df_tcrdock_hcnn = pd.concat([pd.read_csv('hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                 pd.read_csv('hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])

    df_tcrdock_mhc_decoys = pd.read_csv('mhc_pwm/mhc_motif_peptides_w_pae_w_blosum.tsv', sep='\t')

    df_proteinmpnn_002 = pd.read_csv('proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7_w_pae_w_blosum.tsv', sep='\t')
    df_proteinmpnn_02 = pd.read_csv('proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7_w_pae_w_blosum.tsv', sep='\t')

    df_tcrdock_hcnn_fixed_structure = pd.concat([pd.read_csv('hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                                 pd.read_csv('hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])
    
    df_blosum_t1 = pd.read_csv('blosum62/sample_peptides_from_blosum62__temperature=1.0_w_pae_w_blosum.tsv', sep='\t')
    df_blosum_t2 = pd.read_csv('blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae_w_blosum.tsv', sep='\t')
    # df_blosum_t1 = pd.read_csv('proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7_w_pae_w_blosum.tsv', sep='\t')
    # df_blosum_t2 = pd.read_csv('proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7_w_pae_w_blosum.tsv', sep='\t')

    score_columns = ['neg_pmhc_tcr_pae', 'blosum_full', 'blosum_diff', 'neg_pmhc_tcr_pae__blosum_full__min_max__0.5', 'neg_pmhc_tcr_pae__blosum_diff__min_max__0.5']


    ncols = len(score_columns) + 2
    nrows = 1
    colsize = 4
    rowsize = 4
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*colsize, nrows*rowsize), sharex=False, sharey=False)

    for col, score_column in enumerate(score_columns):

        ## histogram of tcrdock neg_pae scores on our designs, colored by HCNN model used to make the designs
        ax = axs[col]
        for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
            scores = df_tcrdock_hcnn.loc[df_tcrdock_hcnn['hcnn_model'] == hcnn_model][score_column].values
            parts = ax.violinplot(scores, positions=[i], widths=0.6, showmedians=True, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(colors[i])
                pc.set_edgecolor(colors[i])
                pc.set_alpha(0.35)
        
        proteinmpnn_002_scores = df_proteinmpnn_002[score_column].values
        parts = ax.violinplot(proteinmpnn_002_scores, positions=[i+1], widths=0.6, showmedians=True, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i+1])
            pc.set_edgecolor(colors[i+1])
            pc.set_alpha(0.35)
        
        proteinmpnn_02_scores = df_proteinmpnn_02[score_column].values
        parts = ax.violinplot(proteinmpnn_02_scores, positions=[i+2], widths=0.6, showmedians=True, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i+2])
            pc.set_edgecolor(colors[i+2])
            pc.set_alpha(0.35)
        
        decoy_scores = df_tcrdock_mhc_decoys[score_column].values
        parts = ax.violinplot(decoy_scores, positions=[i+3], widths=0.6, showmedians=True, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i+3])
            pc.set_edgecolor(colors[i+3])
            pc.set_alpha(0.35)
        
        for j, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
            scores = df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model][score_column].values
            parts = ax.violinplot(scores, positions=[i+4+j], widths=0.6, showmedians=True, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(colors[i+4+j])
                pc.set_edgecolor(colors[i+4+j])
                pc.set_alpha(0.35)
        
        blosum_t1_scores = df_blosum_t1[score_column].values
        parts = ax.violinplot(blosum_t1_scores, positions=[i+4+j+1], widths=0.6, showmedians=True, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i+4+j+1])
            pc.set_edgecolor(colors[i+4+j+1])
            pc.set_alpha(0.35)
        
        blosum_t2_scores = df_blosum_t2[score_column].values
        parts = ax.violinplot(blosum_t2_scores, positions=[i+4+j+2], widths=0.6, showmedians=True, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i+4+j+2])
            pc.set_edgecolor(colors[i+4+j+2])
            pc.set_alpha(0.35)
        
        ax.axhline(df_wildtype.loc[df_wildtype['peptide'] == WT_SEQ][score_column].values[0], color='black', linestyle='--', label=f'WT ({WT_SEQ})')
        ax.set_ylabel(score_column)
        ax.set_xticklabels([])


    ## plot hamming distances
    hcnn_sequences = df_tcrdock_hcnn['peptide']
    hcnn_fixed_structure_sequences = df_tcrdock_hcnn_fixed_structure['peptide']
    hcnn_models = df_tcrdock_hcnn['hcnn_model']

    # plot distribution of hamming distances to wildtype
    hamming_distances = {hcnn_model: [] for hcnn_model in hcnn_models}
    hamming_distances['proteinmpnn 0.02A T=0.7'] = []
    hamming_distances['proteinmpnn 0.2A T=0.7'] = []
    hamming_distances['mhc-sampled decoys'] = []
    for hcnn_model in ['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']:
        hamming_distances[hcnn_model + '_fixed_structure'] = []
    hamming_distances['blosum62 T=1'] = []
    hamming_distances['blosum62 T=2'] = []
    
    for hcnn_model in ['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']:
        sequences = hcnn_sequences[hcnn_models == hcnn_model]
        for seq in sequences:
            hamming_distances[hcnn_model].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
        hamming_distances[hcnn_model] = np.array(hamming_distances[hcnn_model])
    
    proteinmpnn_002_sequences = df_proteinmpnn_002['peptide']
    for seq in proteinmpnn_002_sequences:
        hamming_distances['proteinmpnn 0.02A T=0.7'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
    hamming_distances['proteinmpnn 0.02A T=0.7'] = np.array(hamming_distances['proteinmpnn 0.02A T=0.7'])

    proteinmpnn_02_sequences = df_proteinmpnn_02['peptide']
    for seq in proteinmpnn_02_sequences:
        hamming_distances['proteinmpnn 0.2A T=0.7'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
    hamming_distances['proteinmpnn 0.2A T=0.7'] = np.array(hamming_distances['proteinmpnn 0.2A T=0.7'])
    
    decoy_sequences = df_tcrdock_mhc_decoys['peptide']
    for seq in decoy_sequences:
        hamming_distances['mhc-sampled decoys'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
    hamming_distances['mhc-sampled decoys'] = np.array(hamming_distances['mhc-sampled decoys'])

    for hcnn_model in ['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']:
        sequences = df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model, 'peptide'].values
        for seq in sequences:
            hamming_distances[hcnn_model + '_fixed_structure'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
        hamming_distances[hcnn_model + '_fixed_structure'] = np.array(hamming_distances[hcnn_model + '_fixed_structure'])
    
    for seq in df_blosum_t1['peptide']:
        hamming_distances['blosum62 T=1'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
    hamming_distances['blosum62 T=1'] = np.array(hamming_distances['blosum62 T=1'])

    for seq in df_blosum_t2['peptide']:
        hamming_distances['blosum62 T=2'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
    hamming_distances['blosum62 T=2'] = np.array(hamming_distances['blosum62 T=2'])

    ax = axs[col+1]
    for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        parts = ax.violinplot(hamming_distances[hcnn_model], positions=[i], widths=0.6, showmedians=True, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i])
            pc.set_edgecolor(colors[i])
            pc.set_alpha(0.35)
        ax.text(i, 0.5, str(len(hamming_distances[hcnn_model])), color='black', ha='center', va='center')

    parts = ax.violinplot(hamming_distances['proteinmpnn 0.02A T=0.7'], positions=[i+1], widths=0.6, showmedians=True, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor(colors[i+1])
        pc.set_edgecolor(colors[i+1])
        pc.set_alpha(0.35)
    ax.text(i+1, 0.5, str(len(hamming_distances['proteinmpnn 0.02A T=0.7'])), color='black', ha='center', va='center')

    parts = ax.violinplot(hamming_distances['proteinmpnn 0.2A T=0.7'], positions=[i+2], widths=0.6, showmedians=True, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor(colors[i+2])
        pc.set_edgecolor(colors[i+2])
        pc.set_alpha(0.35)
    ax.text(i+2, 0.5, str(len(hamming_distances['proteinmpnn 0.2A T=0.7'])), color='black', ha='center', va='center')

    parts = ax.violinplot(hamming_distances['mhc-sampled decoys'], positions=[i+3], widths=0.6, showmedians=True, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor(colors[i+3])
        pc.set_edgecolor(colors[i+3])
        pc.set_alpha(0.35)
    ax.text(i+3, 0.5, str(len(hamming_distances['mhc-sampled decoys'])), color='black', ha='center', va='center')

    for j, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        parts = ax.violinplot(hamming_distances[hcnn_model + '_fixed_structure'], positions=[i+4+j], widths=0.6, showmedians=True, showextrema=False)
        for pc in parts['bodies']:
            pc.set_facecolor(colors[i+4+j])
            pc.set_edgecolor(colors[i+4+j])
            pc.set_alpha(0.35)
        ax.text(i+4+j, 0.5, str(len(hamming_distances[hcnn_model + '_fixed_structure'])), color='black', ha='center', va='center')
    
    parts = ax.violinplot(hamming_distances['blosum62 T=1'], positions=[i+4+j+1], widths=0.6, showmedians=True, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor(colors[i+4+j+1])
        pc.set_edgecolor(colors[i+4+j+1])
        pc.set_alpha(0.35)
    ax.text(i+4+j+1, 0.5, str(len(hamming_distances['blosum62 T=1'])), color='black', ha='center', va='center')

    parts = ax.violinplot(hamming_distances['blosum62 T=2'], positions=[i+4+j+2], widths=0.6, showmedians=True, showextrema=False)
    for pc in parts['bodies']:
        pc.set_facecolor(colors[i+4+j+2])
        pc.set_edgecolor(colors[i+4+j+2])
        pc.set_alpha(0.35)
    ax.text(i+4+j+2, 0.5, str(len(hamming_distances['blosum62 T=2'])), color='black', ha='center', va='center')

    ax.set_ylim(0, 9)
    ax.set_ylabel('Hamming distance to wildtype')
    ax.set_xticklabels([])


    ## plot proportion of peptides that are binders according to NetMHCPan
    ax = axs[col+2]
    for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        n_binders = df_tcrdock_hcnn.loc[df_tcrdock_hcnn['hcnn_model'] == hcnn_model]['is_binder_by_netmhc_pan'].sum()
        n_total = len(df_tcrdock_hcnn.loc[df_tcrdock_hcnn['hcnn_model'] == hcnn_model])
        ax.bar(i, n_binders/n_total, color=colors[i], alpha=0.35)
        ax.text(i, 0.05, f'{n_total}', color='black', ha='center', va='center')

    n_binders = df_proteinmpnn_002['is_binder_by_netmhc_pan'].sum()
    n_total = len(df_proteinmpnn_002)
    ax.bar(i+1, n_binders/n_total, color=colors[i+1], alpha=0.35)
    ax.text(i+1, 0.05, f'{n_total}', color='black', ha='center', va='center')

    n_binders = df_proteinmpnn_02['is_binder_by_netmhc_pan'].sum()
    n_total = len(df_proteinmpnn_02)
    ax.bar(i+2, n_binders/n_total, color=colors[i+2], alpha=0.35)
    ax.text(i+2, 0.05, f'{n_total}', color='black', ha='center', va='center')

    n_binders = df_tcrdock_mhc_decoys['is_binder_by_netmhc_pan'].sum()
    n_total = len(df_tcrdock_mhc_decoys)
    ax.bar(i+3, n_binders/n_total, color=colors[i+3], alpha=0.35)
    ax.text(i+3, 0.05, f'{n_total}', color='black', ha='center', va='center')

    for j, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        n_binders = df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model]['is_binder_by_netmhc_pan'].sum()
        n_total = len(df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model])
        ax.bar(i+4+j, n_binders/n_total, color=colors[i+4+j], alpha=0.35)
        ax.text(i+4+j, 0.05, f'{n_total}', color='black', ha='center', va='center')
    
    n_binders = df_blosum_t1['is_binder_by_netmhc_pan'].sum()
    n_total = len(df_blosum_t1)
    ax.bar(i+4+j+1, n_binders/n_total, color=colors[i+4+j+1], alpha=0.35)
    ax.text(i+4+j+1, 0.05, f'{n_total}', color='black', ha='center', va='center')

    n_binders = df_blosum_t2['is_binder_by_netmhc_pan'].sum()
    n_total = len(df_blosum_t2)
    ax.bar(i+4+j+2, n_binders/n_total, color=colors[i+4+j+2], alpha=0.35)
    ax.text(i+4+j+2, 0.05, f'{n_total}', color='black', ha='center', va='center')

    ax.set_ylim(0, 1)
    ax.grid(axis='y', linestyle='--')
    ax.set_ylabel('Proportion of peptides\nthat are binders\naccording to NetMHCPan')
    ax.set_xticklabels([])


    plt.tight_layout()
    plt.savefig('tcrdock_and_blosum_plots.png')
    plt.close()



    ## compute min and max scores for every score column, for plotting purposes
    limits = {}
    for score_col in score_columns:
        limits[score_col] = (np.nanmin([np.nanmin(df_wildtype[score_col]), np.nanmin(df_tcrdock_hcnn[score_col]), np.nanmin(df_tcrdock_mhc_decoys[score_col]), np.nanmin(df_proteinmpnn_002[score_col]), np.nanmin(df_proteinmpnn_02[score_col]), np.nanmin(df_tcrdock_hcnn_fixed_structure[score_col])]),
                            np.nanmax([np.nanmax(df_wildtype[score_col]), np.nanmax(df_tcrdock_hcnn[score_col]), np.nanmax(df_tcrdock_mhc_decoys[score_col]), np.nanmax(df_proteinmpnn_002[score_col]), np.nanmax(df_proteinmpnn_02[score_col]), np.nanmax(df_tcrdock_hcnn_fixed_structure[score_col])]))
        # extend by 5% in each direction
        limits[score_col] = (limits[score_col][0] - 0.05 * (limits[score_col][1] - limits[score_col][0]), limits[score_col][1] + 0.05 * (limits[score_col][1] - limits[score_col][0]))


    ## compute hamming distances

    hcnn_sequences = df_tcrdock_hcnn['peptide']
    hcnn_fixed_structure_sequences = df_tcrdock_hcnn_fixed_structure['peptide']
    hcnn_models = df_tcrdock_hcnn['hcnn_model']

    # plot distribution of hamming distances to wildtype
    hamming_distances = {hcnn_model: [] for hcnn_model in hcnn_models}
    hamming_distances['proteinmpnn 0.02A T=0.7'] = []
    hamming_distances['proteinmpnn 0.2A T=0.7'] = []
    hamming_distances['mhc-sampled decoys'] = []
    for hcnn_model in ['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']:
        hamming_distances[hcnn_model + '_fixed_structure'] = []
    hamming_distances['blosum62 T=1'] = []
    hamming_distances['blosum62 T=2'] = []
    
    for hcnn_model in ['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']:
        sequences = hcnn_sequences[hcnn_models == hcnn_model]
        for seq in sequences:
            hamming_distances[hcnn_model].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
        hamming_distances[hcnn_model] = np.array(hamming_distances[hcnn_model])
    
    proteinmpnn_002_sequences = df_proteinmpnn_002['peptide']
    for seq in proteinmpnn_002_sequences:
        hamming_distances['proteinmpnn 0.02A T=0.7'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
    hamming_distances['proteinmpnn 0.02A T=0.7'] = np.array(hamming_distances['proteinmpnn 0.02A T=0.7'])

    proteinmpnn_02_sequences = df_proteinmpnn_02['peptide']
    for seq in proteinmpnn_02_sequences:
        hamming_distances['proteinmpnn 0.2A T=0.7'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
    hamming_distances['proteinmpnn 0.2A T=0.7'] = np.array(hamming_distances['proteinmpnn 0.2A T=0.7'])
    
    decoy_sequences = df_tcrdock_mhc_decoys['peptide']
    for seq in decoy_sequences:
        hamming_distances['mhc-sampled decoys'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
    hamming_distances['mhc-sampled decoys'] = np.array(hamming_distances['mhc-sampled decoys'])

    for hcnn_model in ['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']:
        sequences = df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model, 'peptide'].values
        for seq in sequences:
            hamming_distances[hcnn_model + '_fixed_structure'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
        hamming_distances[hcnn_model + '_fixed_structure'] = np.array(hamming_distances[hcnn_model + '_fixed_structure'])

    for seq in df_blosum_t1['peptide']:
        hamming_distances['blosum62 T=1'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
    hamming_distances['blosum62 T=1'] = np.array(hamming_distances['blosum62 T=1'])

    for seq in df_blosum_t2['peptide']:
        hamming_distances['blosum62 T=2'].append(sum([seq[i] != WT_SEQ[i] for i in range(len(seq))]))
    hamming_distances['blosum62 T=2'] = np.array(hamming_distances['blosum62 T=2'])


    ## make plots of pae scores, divided by distance to wildtype
    ## histogram of tcrdock neg_pae scores on our designs, colored by HCNN model used to make the designs

    hamming_distances = {}
    hamming_distances['mhc-sampled decoys'] = {}
    for seq in hcnn_sequences:
        hamming_distances[seq] = sum([seq[i] != WT_SEQ[i] for i in range(len(seq))])
    
    for seq in proteinmpnn_002_sequences:
        hamming_distances[seq] = sum([seq[i] != WT_SEQ[i] for i in range(len(seq))])
    
    for seq in proteinmpnn_02_sequences:
        hamming_distances[seq] = sum([seq[i] != WT_SEQ[i] for i in range(len(seq))])
    
    for seq in decoy_sequences:
        hamming_distances[seq] = sum([seq[i] != WT_SEQ[i] for i in range(len(seq))])
    
    for seq in hcnn_fixed_structure_sequences:
        hamming_distances[seq] = sum([seq[i] != WT_SEQ[i] for i in range(len(seq))])
    
    for seq in df_blosum_t1['peptide']:
        hamming_distances[seq] = sum([seq[i] != WT_SEQ[i] for i in range(len(seq))])
    
    for seq in df_blosum_t2['peptide']:
        hamming_distances[seq] = sum([seq[i] != WT_SEQ[i] for i in range(len(seq))])


    ncols = len(WT_SEQ)
    nrows = len(score_columns)
    colsize = 5
    rowsize = 4
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*colsize, nrows*rowsize), sharex=True) #, sharey=True)

    for row, score_col in enumerate(score_columns):

        for col, hamming_distance in enumerate(range(1, len(WT_SEQ)+1)):

            sequences_at_dist = set([seq for seq, dist in hamming_distances.items() if dist == hamming_distance])

            ax = axs[row, col]
            for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
                pep_sequences = df_tcrdock_hcnn.loc[df_tcrdock_hcnn['hcnn_model'] == hcnn_model, 'peptide'].values
                scores = [df_tcrdock_hcnn.loc[df_tcrdock_hcnn['peptide'] == seq][score_col].values[0] for seq in pep_sequences if seq in sequences_at_dist]
                
                if len(scores) == 0:
                    continue

                ## boxplot
                # ax.boxplot(scores, positions=[i], widths=0.6, patch_artist=True, boxprops=dict(facecolor=colors[i], color=colors[i]), medianprops=dict(color='black'), showfliers=False)
                ## violinplot equivalent of the above boxplot command
                parts = ax.violinplot(scores, positions=[i], widths=0.6, showmedians=True, showextrema=False)
                for pc in parts['bodies']:
                    pc.set_facecolor(colors[i])
                    pc.set_edgecolor(colors[i])
                    pc.set_alpha(0.35)
                ## on each violinplot, display the number of points
                ax.text(i, -9.8, str(len(scores)), color='black', ha='center', va='center')

                # ax.hist(scores, bins=20, color=colors[i], alpha=0.35, label=hcnn_model)
                # ax.axvline(np.mean(scores), color=colors[i], linestyle='--')
            
            proteinmpnn_scores = [df_proteinmpnn_002.loc[df_proteinmpnn_002['peptide'] == seq][score_col].values[0] for seq in proteinmpnn_002_sequences if seq in sequences_at_dist]
            if len(proteinmpnn_scores) > 0:
                parts = ax.violinplot(proteinmpnn_scores, positions=[i+1], widths=0.6, showmedians=True, showextrema=False)
                for pc in parts['bodies']:
                    pc.set_facecolor(colors[i+1])
                    pc.set_edgecolor(colors[i+1])
                    pc.set_alpha(0.35)
                ax.text(i+1, -9.8, str(len(proteinmpnn_scores)), color='black', ha='center', va='center')

            proteinmpnn_scores = [df_proteinmpnn_02.loc[df_proteinmpnn_02['peptide'] == seq][score_col].values[0] for seq in proteinmpnn_02_sequences if seq in sequences_at_dist]
            if len(proteinmpnn_scores) > 0:
                parts = ax.violinplot(proteinmpnn_scores, positions=[i+2], widths=0.6, showmedians=True, showextrema=False)
                for pc in parts['bodies']:
                    pc.set_facecolor(colors[i+2])
                    pc.set_edgecolor(colors[i+2])
                    pc.set_alpha(0.35)
                ax.text(i+2, -9.8, str(len(proteinmpnn_scores)), color='black', ha='center', va='center')
            
            decoy_scores = [df_tcrdock_mhc_decoys.loc[df_tcrdock_mhc_decoys['peptide'] == seq][score_col].values[0] for seq in decoy_sequences if seq in sequences_at_dist]
            if len(decoy_scores) > 0:
                parts = ax.violinplot(decoy_scores, positions=[i+3], widths=0.6, showmedians=True, showextrema=False)
                for pc in parts['bodies']:
                    pc.set_facecolor(colors[i+3])
                    pc.set_edgecolor(colors[i+3])
                    pc.set_alpha(0.35)
                ax.text(i+3, -9.8, str(len(decoy_scores)), color='black', ha='center', va='center')

            for k, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
                pep_sequences = df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model, 'peptide'].values
                scores = [df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['peptide'] == seq][score_col].values[0] for seq in pep_sequences if seq in sequences_at_dist]
                if len(scores) == 0:
                    continue

                parts = ax.violinplot(scores, positions=[i+4+k], widths=0.6, showmedians=True, showextrema=False)
                for pc in parts['bodies']:
                    pc.set_facecolor(colors[i+4+k])
                    pc.set_edgecolor(colors[i+4+k])
                    pc.set_alpha(0.35)
                ax.text(i+4+k, -9.8, str(len(scores)), color='black', ha='center', va='center')
            
            blosum_t1_scores = [df_blosum_t1.loc[df_blosum_t1['peptide'] == seq][score_col].values[0] for seq in df_blosum_t1['peptide'] if seq in sequences_at_dist]
            if len(blosum_t1_scores) > 0:
                parts = ax.violinplot(blosum_t1_scores, positions=[i+4+k+1], widths=0.6, showmedians=True, showextrema=False)
                for pc in parts['bodies']:
                    pc.set_facecolor(colors[i+4+k+1])
                    pc.set_edgecolor(colors[i+4+k+1])
                    pc.set_alpha(0.35)
                ax.text(i+4+k+1, -9.8, str(len(blosum_t1_scores)), color='black', ha='center', va='center')
            
            blosum_t2_scores = [df_blosum_t2.loc[df_blosum_t2['peptide'] == seq][score_col].values[0] for seq in df_blosum_t2['peptide'] if seq in sequences_at_dist]
            if len(blosum_t2_scores) > 0:
                parts = ax.violinplot(blosum_t2_scores, positions=[i+4+k+2], widths=0.6, showmedians=True, showextrema=False)
                for pc in parts['bodies']:
                    pc.set_facecolor(colors[i+4+k+2])
                    pc.set_edgecolor(colors[i+4+k+2])
                    pc.set_alpha(0.35)
                ax.text(i+4+k+2, -9.8, str(len(blosum_t2_scores)), color='black', ha='center', va='center')
            
            ax.set_xticklabels([])
            ax.set_title('Hamming distance to WT = {}'.format(hamming_distance))
            ax.set_ylabel(score_col)
            ax.axhline(df_wildtype.loc[df_wildtype['peptide'] == WT_SEQ][score_col].values[0], color='black', linestyle='--', label=f'WT ({WT_SEQ})')
            ax.set_ylim(limits[score_col])
    
    plt.tight_layout()
    plt.savefig('tcrdock_and_blosum_plots_by_hamming_distance_to_WT.png')
    plt.close()



    color_list = plt.get_cmap('tab20').colors
    blue = color_list[0]
    blue_light = color_list[1]
    orange = color_list[2]
    orange_light = color_list[3]
    green = color_list[4]
    green_light = color_list[5]
    red = color_list[6]
    red_light = color_list[7]
    purple = color_list[8]
    purple_light = color_list[9]
    brown = color_list[10]
    brown_light = color_list[11]
    pink = color_list[12]
    pink_light = color_list[13]
    gray = color_list[14]
    gray_light = color_list[15]
    olive = color_list[16]
    olive_light = color_list[17]
    cyan = color_list[18]
    cyan_light = color_list[19]

    colors_here = [red, red_light, orange, orange_light, green, green_light, brown]

    hamming_distances_to_consider = [3, 4, 5, 6, 7, 8, 9]

    text_y = -8.6

    score_col = 'neg_pmhc_tcr_pae'
    ncols = len(hamming_distances_to_consider)
    nrows = 1
    colsize = 5
    rowsize = 4
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*colsize, nrows*rowsize), sharex=True, sharey=True)

    for col, hamming_distance in enumerate(hamming_distances_to_consider):

        sequences_at_dist = set([seq for seq, dist in hamming_distances.items() if dist == hamming_distance])

        ax = axs[col]

        for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
            pep_sequences = df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model, 'peptide'].values
            scores = [df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['peptide'] == seq][score_col].values[0] for seq in pep_sequences if seq in sequences_at_dist]
            if len(scores) == 0:
                continue

            parts = ax.violinplot(scores, positions=[i], widths=0.6, showmedians=False, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(colors_here[i])
                pc.set_edgecolor(colors_here[i])
                pc.set_alpha(0.5)
            median = np.median(scores)
            ax.plot([i - 0.15, i + 0.15], [median, median], color=colors_here[i], lw=2)
            ax.text(i, text_y, str(len(scores)), color='black', ha='center', va='center')

        for k, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
            pep_sequences = df_tcrdock_hcnn.loc[df_tcrdock_hcnn['hcnn_model'] == hcnn_model, 'peptide'].values
            scores = [df_tcrdock_hcnn.loc[df_tcrdock_hcnn['peptide'] == seq][score_col].values[0] for seq in pep_sequences if seq in sequences_at_dist]
            
            if len(scores) == 0:
                continue

            ## boxplot
            # ax.boxplot(scores, positions=[i], widths=0.6, patch_artist=True, boxprops=dict(facecolor=colors[i], color=colors[i]), medianprops=dict(color='black'), showfliers=False)
            ## violinplot equivalent of the above boxplot command
            parts = ax.violinplot(scores, positions=[i+k+1], widths=0.6, showmedians=False, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(colors_here[i+k+1])
                pc.set_edgecolor(colors_here[i+k+1])
                pc.set_alpha(0.5)
            median = np.median(scores)
            ax.plot([i+k+1 - 0.15, i+k+1 + 0.15], [median, median], color=colors_here[i+k+1], lw=2)
            ## on each violinplot, display the number of points
            ax.text(i+k+1, text_y, str(len(scores)), color='black', ha='center', va='center')

            # ax.hist(scores, bins=20, color=colors_here[i], alpha=0.35, label=hcnn_model)
            # ax.axvline(np.mean(scores), color=colors_here[i], linestyle='--')
        
        proteinmpnn_scores = [df_proteinmpnn_002.loc[df_proteinmpnn_002['peptide'] == seq][score_col].values[0] for seq in proteinmpnn_002_sequences if seq in sequences_at_dist]
        if len(proteinmpnn_scores) > 0:
            parts = ax.violinplot(proteinmpnn_scores, positions=[i+k+2], widths=0.6, showmedians=False, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(colors_here[i+k+2])
                pc.set_edgecolor(colors_here[i+k+2])
                pc.set_alpha(0.5)
            median = np.median(proteinmpnn_scores)
            ax.plot([i+k+2 - 0.15, i+k+2 + 0.15], [median, median], color=colors_here[i+k+2], lw=2)
            ax.text(i+k+2, text_y, str(len(proteinmpnn_scores)), color='black', ha='center', va='center')

        proteinmpnn_scores = [df_proteinmpnn_02.loc[df_proteinmpnn_02['peptide'] == seq][score_col].values[0] for seq in proteinmpnn_02_sequences if seq in sequences_at_dist]
        if len(proteinmpnn_scores) > 0:
            parts = ax.violinplot(proteinmpnn_scores, positions=[i+k+3], widths=0.6, showmedians=False, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(colors_here[i+k+3])
                pc.set_edgecolor(colors_here[i+k+3])
                pc.set_alpha(0.5)
            median = np.median(proteinmpnn_scores)
            ax.plot([i+k+3 - 0.15, i+k+3 + 0.15], [median, median], color=colors_here[i+k+3], lw=2)
            ax.text(i+k+3, text_y, str(len(proteinmpnn_scores)), color='black', ha='center', va='center')
        
        decoy_scores = [df_tcrdock_mhc_decoys.loc[df_tcrdock_mhc_decoys['peptide'] == seq][score_col].values[0] for seq in decoy_sequences if seq in sequences_at_dist]
        if len(decoy_scores) > 0:
            parts = ax.violinplot(decoy_scores, positions=[i+k+4], widths=0.6, showmedians=False, showextrema=False)
            for pc in parts['bodies']:
                pc.set_facecolor(colors_here[i+k+4])
                pc.set_edgecolor(colors_here[i+k+4])
                pc.set_alpha(0.5)
            median = np.median(decoy_scores)
            ax.plot([i+k+4 - 0.15, i+k+4 + 0.15], [median, median], color=colors_here[i+k+4], lw=2)
            ax.text(i+k+4, text_y, str(len(decoy_scores)), color='black', ha='center', va='center')
        
        
        ax.set_xticklabels([])
        ax.set_title('Hamming distance to WT = {}'.format(hamming_distance))
        if col == 0: ax.set_ylabel(score_col)
        ax.axhline(df_wildtype.loc[df_wildtype['peptide'] == WT_SEQ][score_col].values[0], color='black', linestyle='--', label=f'WT ({WT_SEQ})')
        ax.set_ylim(limits[score_col])
    
    plt.tight_layout()
    plt.savefig('__pretty_tcrdock_and_blosum_plots_by_hamming_distance_to_WT.png')
    plt.savefig('__pretty_tcrdock_and_blosum_plots_by_hamming_distance_to_WT.pdf')
    plt.close()


    score_col = 'is_binder_by_netmhc_pan'
    plt.figure(figsize=(4, 3))

    bars = []

    for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        scores = df_tcrdock_hcnn_fixed_structure[score_col].values
        bars.append(np.mean(scores) * 100)

    for k, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        scores = df_tcrdock_hcnn[score_col].values
        bars.append(np.mean(scores) * 100)
    
    proteinmpnn_scores = df_proteinmpnn_002[score_col].values
    bars.append(np.mean(proteinmpnn_scores) * 100)

    proteinmpnn_scores = df_proteinmpnn_02[score_col].values
    bars.append(np.mean(proteinmpnn_scores) * 100)

    decoy_scores = df_tcrdock_mhc_decoys[score_col].values
    bars.append(np.mean(decoy_scores) * 100)

    plt.bar(range(len(bars)), bars, color=colors_here, alpha=0.8)
    plt.axhline(100, color='black', linestyle='-')

    plt.ylabel('Percent of peptides\nthat bind MHC (NetMHCPan)')
    
    plt.tight_layout()
    plt.savefig('__pretty_is_netmhc_binder.png')
    plt.savefig('__pretty_is_netmhc_binder.pdf')
    plt.close()


    score_col = 'is_self_peptide'
    plt.figure(figsize=(4, 3))

    bars = []

    for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        scores = df_tcrdock_hcnn_fixed_structure[score_col].values
        bars.append(np.mean(scores) * 100)

    for k, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        scores = df_tcrdock_hcnn[score_col].values
        bars.append(np.mean(scores) * 100)
    
    proteinmpnn_scores = df_proteinmpnn_002[score_col].values
    bars.append(np.mean(proteinmpnn_scores) * 100)

    proteinmpnn_scores = df_proteinmpnn_02[score_col].values
    bars.append(np.mean(proteinmpnn_scores) * 100)

    decoy_scores = df_tcrdock_mhc_decoys[score_col].values
    bars.append(np.mean(decoy_scores) * 100)

    plt.bar(range(len(bars)), bars, color=colors_here, alpha=0.8)
    plt.axhline(100, color='black', linestyle='-')

    plt.ylabel('Percent of peptides\nthat are self peptides.')
    
    plt.tight_layout()
    plt.savefig('__pretty_is_self.png')
    plt.savefig('__pretty_is_self.pdf')
    plt.close()



    fig, ax = plt.subplots(figsize=(4, 4))
    handles = []
    for j, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        handles.append(Patch(facecolor=colors_here[j], edgecolor=colors_here[j], alpha=0.5, label=hcnn_model + '__fixed_structure'))
    for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        handles.append(Patch(facecolor=colors_here[i+j], edgecolor=colors_here[i+j], alpha=0.5, label=hcnn_model + '__sim_anneal'))
    handles.append(Patch(facecolor=colors_here[i+j+1], edgecolor=colors_here[i+j+1], alpha=0.5, label='proteinmpnn 0.02A T=0.7'))
    handles.append(Patch(facecolor=colors_here[i+j+2], edgecolor=colors_here[i+j+2], alpha=0.5, label='proteinmpnn 0.2A T=0.7'))
    handles.append(Patch(facecolor=colors_here[i+j+3], edgecolor=colors_here[i+j+3], alpha=0.5, label='mhc-sampled decoys'))
    handles.append(Line2D([0], [0], color='black', linestyle='--', label=f'WT ({WT_SEQ})'))
    ax.legend(handles=handles, loc='center')
    plt.axis('off')
    plt.savefig('__pretty_tcrdock_and_blosum_plots_legend_patches.png')
    plt.savefig('__pretty_tcrdock_and_blosum_plots_legend_patches.pdf')
    plt.close()



    fig, ax = plt.subplots(figsize=(4, 4))
    handles = []
    for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        handles.append(Patch(facecolor=colors[i], edgecolor=colors[i], alpha=0.35, label=hcnn_model))
    handles.append(Patch(facecolor=colors[i+1], edgecolor=colors[i+1], alpha=0.35, label='proteinmpnn 0.02A T=0.7'))
    handles.append(Patch(facecolor=colors[i+2], edgecolor=colors[i+2], alpha=0.35, label='proteinmpnn 0.2A T=0.7'))
    handles.append(Patch(facecolor=colors[i+3], edgecolor=colors[i+3], alpha=0.35, label='mhc-sampled decoys'))
    for j, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        handles.append(Patch(facecolor=colors[i+4+j], edgecolor=colors[i+4+j], alpha=0.35, label=hcnn_model + '_fixed_structure'))
    handles.append(Patch(facecolor=colors[i+4+j+1], edgecolor=colors[i+4+j+1], alpha=0.35, label='blosum62 T=1'))
    handles.append(Patch(facecolor=colors[i+4+j+2], edgecolor=colors[i+4+j+2], alpha=0.35, label='blosum62 T=2'))
    handles.append(Line2D([0], [0], color='black', linestyle='--', label=f'WT ({WT_SEQ})'))
    ax.legend(handles=handles, loc='center')
    plt.axis('off')
    plt.savefig('tcrdock_and_blosum_plots_legend_patches.png')
    plt.close()


    # make logoplot of highest-scoring peptides
    import logomaker

    cutoff = -5.5

    nrows = 11
    ncols = 1
    rowsize = 1.5
    colsize = 8
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*colsize, nrows*rowsize), sharex=True, sharey=True)
    
    def make_pwm_df(sequences):
        pwm = np.zeros((len(WT_SEQ), 20))
        for seq in sequences:
            assert len(seq) == len(WT_SEQ)
            for i, aa in enumerate(seq):
                pwm[i, [ol_to_ind_size[aa]]] += 1
        for j in range(len(WT_SEQ)):
            pwm[j, :] /= pwm[j, :].sum()
        
        aas, indices = zip(*sorted(list(ol_to_ind_size.items()), key=lambda x: x[1]))
        df_logo = pd.DataFrame(pwm, columns=aas)
        # make nan values zero
        df_logo.fillna(0, inplace=True)
        return df_logo

    for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        scores = -df_tcrdock_hcnn.loc[df_tcrdock_hcnn['hcnn_model'] == hcnn_model]['pmhc_tcr_pae'].values
        good_enough_seqs = df_tcrdock_hcnn.loc[df_tcrdock_hcnn['hcnn_model'] == hcnn_model]['peptide'].values[scores > cutoff]

        df_logo = make_pwm_df(good_enough_seqs)
        ax = axs[i]
        logomaker.Logo(df_logo, ax=ax)
        ax.set_title(hcnn_model + f' - {len(good_enough_seqs)} peptides below {-cutoff} PAE', fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])

    proteinmpnn_scores = -df_proteinmpnn_002['pmhc_tcr_pae'].values
    good_enough_seqs = df_proteinmpnn_002['peptide'].values[proteinmpnn_scores > cutoff]
    df_logo = make_pwm_df(good_enough_seqs)
    ax = axs[i+1]
    logomaker.Logo(df_logo, ax=ax)
    ax.set_title('proteinmpnn 0.02A T=0.7' + f' - {len(good_enough_seqs)} peptides below {-cutoff} PAE', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    
    proteinmpnn_scores = -df_proteinmpnn_02['pmhc_tcr_pae'].values
    good_enough_seqs = df_proteinmpnn_02['peptide'].values[proteinmpnn_scores > cutoff]
    df_logo = make_pwm_df(good_enough_seqs)
    ax = axs[i+2]
    logomaker.Logo(df_logo, ax=ax)
    ax.set_title('proteinmpnn 0.2A T=0.7' + f' - {len(good_enough_seqs)} peptides below {-cutoff} PAE', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])

    decoy_scores = -df_tcrdock_mhc_decoys['pmhc_tcr_pae'].values
    good_enough_seqs = df_tcrdock_mhc_decoys['peptide'].values[decoy_scores > cutoff]
    df_logo = make_pwm_df(good_enough_seqs)
    ax = axs[i+3]
    logomaker.Logo(df_logo, ax=ax)
    ax.set_title('mhc-sampled decoys' + f' - {len(good_enough_seqs)} peptides below {-cutoff} PAE', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])

    for j, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        scores = -df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model]['pmhc_tcr_pae'].values
        good_enough_seqs = df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model]['peptide'].values[scores > cutoff]
        df_logo = make_pwm_df(good_enough_seqs)
        ax = axs[i+4+j]
        logomaker.Logo(df_logo, ax=ax)
        ax.set_title(hcnn_model + '_fixed_structure' + f' - {len(good_enough_seqs)} peptides below {-cutoff} PAE', fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])
    
    blosum_t1_scores = -df_blosum_t1['pmhc_tcr_pae'].values
    good_enough_seqs = df_blosum_t1['peptide'].values[blosum_t1_scores > cutoff]
    df_logo = make_pwm_df(good_enough_seqs)
    ax = axs[i+4+j+1]
    logomaker.Logo(df_logo, ax=ax)
    ax.set_title('blosum62 T=1' + f' - {len(good_enough_seqs)} peptides below {-cutoff} PAE', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])

    blosum_t2_scores = -df_blosum_t2['pmhc_tcr_pae'].values
    good_enough_seqs = df_blosum_t2['peptide'].values[blosum_t2_scores > cutoff]
    df_logo = make_pwm_df(good_enough_seqs)
    ax = axs[i+4+j+2]
    logomaker.Logo(df_logo, ax=ax)
    ax.set_title('blosum62 T=2' + f' - {len(good_enough_seqs)} peptides below {-cutoff} PAE', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(f'tcrdock_logos_of_pae_scores_below_{-cutoff}.png')
    plt.close()

    

    ## same as above, but for bad peptides

    cutoff = -7

    nrows = 11
    ncols = 1
    rowsize = 1.5
    colsize = 8
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(ncols*colsize, nrows*rowsize), sharex=True, sharey=True)

    for i, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        scores = -df_tcrdock_hcnn.loc[df_tcrdock_hcnn['hcnn_model'] == hcnn_model]['pmhc_tcr_pae'].values
        bad_enough_seqs = df_tcrdock_hcnn.loc[df_tcrdock_hcnn['hcnn_model'] == hcnn_model]['peptide'].values[scores < cutoff]

        df_logo = make_pwm_df(bad_enough_seqs)
        ax = axs[i]
        logomaker.Logo(df_logo, ax=ax)
        ax.set_title(hcnn_model + f' - {len(bad_enough_seqs)} peptides above {-cutoff} PAE', fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])
    
    proteinmpnn_scores = -df_proteinmpnn_002['pmhc_tcr_pae'].values
    bad_enough_seqs = df_proteinmpnn_002['peptide'].values[proteinmpnn_scores < cutoff]
    df_logo = make_pwm_df(bad_enough_seqs)
    ax = axs[i+1]
    logomaker.Logo(df_logo, ax=ax)
    ax.set_title('proteinmpnn 0.02A T=0.7' + f' - {len(bad_enough_seqs)} peptides above {-cutoff} PAE', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])
    
    proteinmpnn_scores = -df_proteinmpnn_02['pmhc_tcr_pae'].values
    bad_enough_seqs = df_proteinmpnn_02['peptide'].values[proteinmpnn_scores < cutoff]
    df_logo = make_pwm_df(bad_enough_seqs)
    ax = axs[i+2]
    logomaker.Logo(df_logo, ax=ax)
    ax.set_title('proteinmpnn 0.2A T=0.7' + f' - {len(bad_enough_seqs)} peptides above {-cutoff} PAE', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])

    decoy_scores = -df_tcrdock_mhc_decoys['pmhc_tcr_pae'].values
    bad_enough_seqs = df_tcrdock_mhc_decoys['peptide'].values[decoy_scores < cutoff]
    df_logo = make_pwm_df(bad_enough_seqs)
    ax = axs[i+3]
    logomaker.Logo(df_logo, ax=ax)
    ax.set_title('mhc-sampled decoys' + f' - {len(bad_enough_seqs)} peptides above {-cutoff} PAE', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])

    for j, hcnn_model in enumerate(['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']):
        scores = -df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model]['pmhc_tcr_pae'].values
        bad_enough_seqs = df_tcrdock_hcnn_fixed_structure.loc[df_tcrdock_hcnn_fixed_structure['hcnn_model'] == hcnn_model]['peptide'].values[scores < cutoff]
        df_logo = make_pwm_df(bad_enough_seqs)
        ax = axs[i+4+j]
        logomaker.Logo(df_logo, ax=ax)
        ax.set_title(hcnn_model + '_fixed_structure' + f' - {len(bad_enough_seqs)} peptides above {-cutoff} PAE', fontsize=14)
        ax.set_xticks([])
        ax.set_yticks([])

    blosum_t1_scores = -df_blosum_t1['pmhc_tcr_pae'].values
    bad_enough_seqs = df_blosum_t1['peptide'].values[blosum_t1_scores < cutoff]
    df_logo = make_pwm_df(bad_enough_seqs)
    ax = axs[i+4+j+1]
    logomaker.Logo(df_logo, ax=ax)
    ax.set_title('blosum62 T=1' + f' - {len(bad_enough_seqs)} peptides above {-cutoff} PAE', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])

    blosum_t2_scores = -df_blosum_t2['pmhc_tcr_pae'].values
    bad_enough_seqs = df_blosum_t2['peptide'].values[blosum_t2_scores < cutoff]
    df_logo = make_pwm_df(bad_enough_seqs)
    ax = axs[i+4+j+2]
    logomaker.Logo(df_logo, ax=ax)
    ax.set_title('blosum62 T=2' + f' - {len(bad_enough_seqs)} peptides above {-cutoff} PAE', fontsize=14)
    ax.set_xticks([])
    ax.set_yticks([])

    plt.tight_layout()
    plt.savefig(f'tcrdock_logos_of_pae_scores_above_{-cutoff}.png')
    plt.close()

            

