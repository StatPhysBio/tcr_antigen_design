
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams['figure.facecolor'] = 'white'

from scipy.stats import spearmanr

THIS_DIR = os.path.dirname(os.path.abspath(__file__))

SYSTEMS = ['nyeso', 'tax', 'hsiue', 'tcr1', 'tcr2', 'tcr3', 'tcr4', 'tcr5', 'tcr6', 'tcr7']

PROTEINMPNN_MODEL_GROUPS = {'proteinmpnn_v_48_002': ['proteinmpnn_v_48_002', 'proteinmpnn_v_48_002_full_pep_masked', 'proteinmpnn_v_48_002_full_pep_and_tcr_masked', 'proteinmpnn_v_48_002_full_tcr_masked'],
                            'proteinmpnn_v_48_020': ['proteinmpnn_v_48_020', 'proteinmpnn_v_48_020_full_pep_masked', 'proteinmpnn_v_48_020_full_pep_and_tcr_masked', 'proteinmpnn_v_48_020_full_tcr_masked']}

def make_pretty_name(model):
    if model == 'proteinmpnn_v_48_002':
        return 'Pep 1-by-1 Masked'
    elif model == 'proteinmpnn_v_48_002_full_pep_masked':
        return 'Pep Full Masked'
    elif model == 'proteinmpnn_v_48_002_full_tcr_masked':
        return 'TCR Full Masked'
    elif model == 'proteinmpnn_v_48_002_full_pep_and_tcr_masked':
        return 'Pep and TCR Full Masked'
    elif model == 'proteinmpnn_v_48_020':
        return 'Pep 1-by-1 Masked'
    elif model == 'proteinmpnn_v_48_020_full_pep_masked':
        return 'Pep Full Masked'
    elif model == 'proteinmpnn_v_48_020_full_tcr_masked':
        return 'TCR Full Masked'
    elif model == 'proteinmpnn_v_48_020_full_pep_and_tcr_masked':
        return 'Pep and TCR Full Masked'
    else:
        raise ValueError(f'Unknown model: {model}')

PRED_COL = 'pnlogp'

def system_to_dir(system):

    if system == 'nyeso':
        systemdir = os.path.join(THIS_DIR, '..', 'nyeso')
    elif system == 'tax':
        systemdir = os.path.join(THIS_DIR, '..', 'tax')
    elif system == 'hsiue':
        systemdir = os.path.join(THIS_DIR, '..', 'hsiue_et_al')
    elif 'tcr' in system:
        systemdir = os.path.join(THIS_DIR, '..', 'mskcc')
    
    return systemdir

def system_to_predictions_df(system, proteinmpnn_model_and_protocol):

    if system == 'nyeso':
        df = pd.read_csv(os.path.join(THIS_DIR, '..', 'nyeso', 'results', proteinmpnn_model_and_protocol, 'zero_shot_predictions', 'nyeso_peptide_kd_closest-num_seq_per_target=10-use_mt_structure=0.csv'))
    elif system == 'tax':
        df = pd.read_csv(os.path.join(THIS_DIR, '..', 'tax', 'results', proteinmpnn_model_and_protocol, 'zero_shot_predictions', 'tax_peptide_kd_closest-num_seq_per_target=10-use_mt_structure=0.csv'))
    elif system == 'hsiue':
        df = pd.read_csv(os.path.join(THIS_DIR, '..', 'hsiue_et_al', 'results', proteinmpnn_model_and_protocol, 'zero_shot_predictions', 'hsiue_et_al_H2_sat_mut-num_seq_per_target=10-use_mt_structure=0.csv'))
        
        # only keep one wildtype row
        is_wt_mask = df['is_wt'].astype(bool)
        first_wt_mask = np.zeros(df.shape[0], dtype=bool)
        first_wt_mask[np.arange(df.shape[0])[is_wt_mask][0]] = True
        mask = np.logical_or(~is_wt_mask, first_wt_mask)

        df = df[mask]

    elif 'tcr' in system:
        df = pd.read_csv(os.path.join(THIS_DIR, '..', 'mskcc', 'results', proteinmpnn_model_and_protocol, 'zero_shot_predictions', f'mskcc_{system}_ec50_sat_mut_af3-num_seq_per_target=10-use_mt_structure=0.csv'))

        # only keep one wildtype row
        is_wt_mask = df['is_wt'].astype(bool)
        first_wt_mask = np.zeros(df.shape[0], dtype=bool)
        first_wt_mask[np.arange(df.shape[0])[is_wt_mask][0]] = True
        mask = np.logical_or(~is_wt_mask, first_wt_mask)

        # only keep reliable rows
        mask = np.logical_and(df['is_reliable'].values, mask)

        df = df[mask]

    return df


for system in SYSTEMS:

    systemdir = system_to_dir(system)

    for proteinmpnn_model_version in PROTEINMPNN_MODEL_GROUPS:

        models_list = PROTEINMPNN_MODEL_GROUPS[proteinmpnn_model_version]

        model_to_predictions = {}
        for proteinmpnn_model_and_protocol in models_list:
            df = system_to_predictions_df(system, proteinmpnn_model_and_protocol)
            predictions = df[PRED_COL]
            model_to_predictions[proteinmpnn_model_and_protocol] = predictions
        
        N = len(models_list)
        corr_heatmap = np.full((N, N), np.nan)
        pval_heatmap = np.full((N, N), np.nan)

        for i in range(N):
            for j in range(i):
                model1 = models_list[i]
                model2 = models_list[j]
                sr, sr_pval = spearmanr(model_to_predictions[model1], model_to_predictions[model2])
                corr_heatmap[i, j] = sr
                pval_heatmap[i, j] = sr_pval
        
        # only keep the bottom-left (N-1, N-1) matrices
        corr_heatmap = corr_heatmap[1:N, :N-1]
        pval_heatmap = pval_heatmap[1:N, :N-1]

        df_heatmap = pd.DataFrame(corr_heatmap)

        xticklabels = [make_pretty_name(model) for model in models_list[:-1]]
        yticklabels = [make_pretty_name(model) for model in models_list[1:]]

        if proteinmpnn_model_version == 'proteinmpnn_v_48_002':
            title = 'ProteinMPNN 0.02'
        elif proteinmpnn_model_version == 'proteinmpnn_v_48_020':
            title = 'ProteinMPNN 0.20'
        else:
            raise ValueError(f'Unknown model version: {proteinmpnn_model_version}')

        fontsize = 14

        plt.figure(figsize=(5.5, 5.5))
        ax = sns.heatmap(df_heatmap, cmap='coolwarm', vmin=-1, vmax=1, annot=False)

        # Annotate each cell manually
        for i in range(df_heatmap.shape[0]):
            for j in range(i+1):
                value = df_heatmap.iloc[i, j]
                if pval_heatmap[i, j] >= 0.05:
                    color = 'red'
                else:
                    if np.abs(value) < 0.5:
                        color = 'black'
                    else:
                        color = 'white'
                ax.text(j + 0.5, i + 0.5, f"{value:.2f}", ha='center', va='center', color=color, fontsize=fontsize)

        plt.xticks(np.arange(len(xticklabels))+0.5, xticklabels, rotation=70, fontsize=fontsize, ha='right')
        plt.yticks(np.arange(len(yticklabels))+0.5, yticklabels, rotation=0, fontsize=fontsize)
        plt.title(title, fontsize=fontsize)
        plt.tight_layout()
        plt.savefig(os.path.join(systemdir, f'{system}_{proteinmpnn_model_version}__scoring_protocol_spearman_correlations.png'))
        plt.savefig(os.path.join(systemdir, f'{system}_{proteinmpnn_model_version}__scoring_protocol_spearman_correlations.pdf'))
        plt.close()
        
