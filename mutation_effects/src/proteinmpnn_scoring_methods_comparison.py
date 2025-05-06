
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

SYSTEMS = ['nyeso', 'tax', 'tcr1', 'tcr2', 'tcr3', 'tcr4', 'tcr5', 'tcr6', 'tcr7', 'hsiue']

SYSTEM_TO_COLOR = {
    'nyeso': 'tab:blue',
    'tax': 'tab:orange',
    'hsiue': 'tab:green',
    'tcr1': 'tab:red',
    'tcr2': 'tab:purple',
    'tcr3': 'tab:brown',
    'tcr4': 'tab:pink',
    'tcr5': 'tab:gray',
    'tcr6': 'tab:olive',
    'tcr7': 'tab:cyan',
}

SYSTEM_TO_MARKER = {
    'nyeso': 'o',
    'tax': 's',
    'hsiue': 'd',
    'tcr1': '*',
    'tcr2': '^',
    'tcr3': 'v',
    'tcr4': '<',
    'tcr5': '>',
    'tcr6': 'p',
    'tcr7': '8',
}

SYSTEM_TO_PRETTY_NAME = {
    'nyeso': '1G4 TCR',
    'tax': 'A6 TCR',
    'hsiue': 'H2-scDb (Fab)',
    'tcr1': 'TCR1',
    'tcr2': 'TCR2',
    'tcr3': 'TCR3',
    'tcr4': 'TCR4',
    'tcr5': 'TCR5',
    'tcr6': 'TCR6',
    'tcr7': 'TCR7'
}

PROTEINMPNN_MODEL_GROUPS = {'proteinmpnn_v_48_002': ['proteinmpnn_v_48_002', 'proteinmpnn_v_48_002_full_pep_masked', 'proteinmpnn_v_48_002_full_pep_and_tcr_masked', 'proteinmpnn_v_48_002_full_tcr_masked'],
                            'proteinmpnn_v_48_020': ['proteinmpnn_v_48_020', 'proteinmpnn_v_48_020_full_pep_masked', 'proteinmpnn_v_48_020_full_pep_and_tcr_masked', 'proteinmpnn_v_48_020_full_tcr_masked']}

def make_pretty_name(model):
    if model == 'proteinmpnn_v_48_002':
        return 'Pep 1-by-1'
    elif model == 'proteinmpnn_v_48_002_full_pep_masked':
        return 'Pep'
    elif model == 'proteinmpnn_v_48_002_full_tcr_masked':
        return 'TCR'
    elif model == 'proteinmpnn_v_48_002_full_pep_and_tcr_masked':
        return 'Pep and TCR'
    elif model == 'proteinmpnn_v_48_020':
        return 'Pep 1-by-1'
    elif model == 'proteinmpnn_v_48_020_full_pep_masked':
        return 'Pep'
    elif model == 'proteinmpnn_v_48_020_full_tcr_masked':
        return 'TCR'
    elif model == 'proteinmpnn_v_48_020_full_pep_and_tcr_masked':
        return 'Pep and TCR'
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

def system_to_target_col(system):

    if system == 'nyeso':
        target_col = '-log10(Kd)'
    elif system == 'tax':
        target_col = '-log10(Kd)'
    elif system == 'hsiue':
        target_col = 'IFN_gamma (pg/ml)'
    elif 'tcr' in system:
        target_col = '-log10(EC50)'
    
    return target_col


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

all_correlations = {}

for system in SYSTEMS:

    all_correlations[system] = {}

    systemdir = system_to_dir(system)

    for proteinmpnn_model_version in PROTEINMPNN_MODEL_GROUPS:

        all_correlations[system][proteinmpnn_model_version] = {}

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

                all_correlations[system][proteinmpnn_model_version][(model1, model2)] = (sr, sr_pval)
        
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



## try to put all these correlations in one scatterplot,because the heatmaps are a little bit too much.
## on one axis, we have the 6 possible scoring protocols combinations
## on the other axis, spearman correlation between the two scoring protocols
## each dot is a dataset. Just color each dataset differently, make a legend for it

## Also, make a similar plot, with the same colors, where on the first axis there are the three scoring protocols, vs. the ground truth
## put the two plots next to each other

## make the dots more transparent if they aren't significant (p-value > 0.05)

from itertools import combinations

for proteinmpnn_model_version in PROTEINMPNN_MODEL_GROUPS:

    protocols = PROTEINMPNN_MODEL_GROUPS[proteinmpnn_model_version]

    protocol_pairs = list(combinations(protocols, 2))

    plt.figure(figsize=(4.5, 5))
    ax = plt.gca()

    for i, protocol_pair in enumerate(protocol_pairs):

        for system in SYSTEMS:

            try:
                sr, sr_pval = all_correlations[system][proteinmpnn_model_version][protocol_pair]
            except:
                sr, sr_pval = all_correlations[system][proteinmpnn_model_version][(protocol_pair[1], protocol_pair[0])]

            alpha = 0.95 if sr_pval < 0.05 else 0.2
            color = SYSTEM_TO_COLOR[system]
            marker = SYSTEM_TO_MARKER[system]

            ax.scatter(i, sr, color=color, marker=marker, alpha=alpha, s=90)

    protocol_pairs_pretty_names = [f'{make_pretty_name(model1)} vs. {make_pretty_name(model2)}'  for model1, model2 in protocol_pairs]

    ax.grid(axis='y', ls='--', alpha=0.5)
    ax.set_xticks(np.arange(len(protocol_pairs_pretty_names)))
    ax.set_xticklabels(protocol_pairs_pretty_names, ha='right')
    ax.tick_params(axis='x', labelsize=fontsize-1, rotation=70)
    ax.tick_params(axis='y', labelsize=fontsize-1)

    ax.set_ylabel('Spearman r', fontsize=fontsize)

    if proteinmpnn_model_version == 'proteinmpnn_v_48_002':
        title = 'ProteinMPNN 0.02'
    elif proteinmpnn_model_version == 'proteinmpnn_v_48_020':
        title = 'ProteinMPNN 0.20'
    else:
        raise ValueError(f'Unknown model version: {proteinmpnn_model_version}')

    ax.set_title(title, fontsize=fontsize+1)

    ax.set_ylim((-1.05, 1.05))
    ax.axhline(0, ls='--', color='black')

    plt.tight_layout()
    plt.savefig(f'../scoring_scheme_comparison__{proteinmpnn_model_version}.png')
    plt.savefig(f'../scoring_scheme_comparison__{proteinmpnn_model_version}.pdf')
    plt.close()


for proteinmpnn_model_version in PROTEINMPNN_MODEL_GROUPS:

    protocols = PROTEINMPNN_MODEL_GROUPS[proteinmpnn_model_version]

    plt.figure(figsize=(4, 4.4))
    ax = plt.gca()

    for i, protocol in enumerate(protocols):

        for system in SYSTEMS:

            df = system_to_predictions_df(system, protocol)

            sr, sr_pval = spearmanr(df[PRED_COL], df[system_to_target_col(system)])


            alpha = 0.95 if sr_pval < 0.05 else 0.2
            color = SYSTEM_TO_COLOR[system]
            marker = SYSTEM_TO_MARKER[system]

            ax.scatter(i, sr, color=color, marker=marker, alpha=alpha, s=90)

    protocol_pretty_names = [make_pretty_name(model) + ' vs. Exp' for model in protocols]

    ax.set_xticks(np.arange(len(protocol_pretty_names)))
    ax.set_xticklabels(protocol_pretty_names, ha='right')
    ax.tick_params(axis='x', labelsize=fontsize-1, rotation=70)
    ax.tick_params(axis='y', labelsize=fontsize-1)

    ax.grid(axis='y', ls='--', alpha=0.5)

    ax.set_ylabel('Spearman r', fontsize=fontsize)

    if proteinmpnn_model_version == 'proteinmpnn_v_48_002':
        title = 'ProteinMPNN 0.02'
    elif proteinmpnn_model_version == 'proteinmpnn_v_48_020':
        title = 'ProteinMPNN 0.20'
    else:
        raise ValueError(f'Unknown model version: {proteinmpnn_model_version}')

    ax.set_title(title, fontsize=fontsize+1)

    ax.set_ylim((-1.05, 1.05))
    ax.axhline(0, ls='--', color='black')

    plt.tight_layout()
    plt.savefig(f'../scoring_schemes_vs_true__{proteinmpnn_model_version}.png')
    plt.savefig(f'../scoring_schemes_vs_true__{proteinmpnn_model_version}.pdf')
    plt.close()




# make legend
# Create dummy legend handles
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
legend_handles = []
for system in SYSTEMS:
    handle = mlines.Line2D([], [], 
                           color=SYSTEM_TO_COLOR[system], 
                           marker=SYSTEM_TO_MARKER[system], 
                           linestyle='None',
                           markersize=10,
                           label=SYSTEM_TO_PRETTY_NAME[system])
    legend_handles.append(handle)

# Create a legend-only plot
fig, ax = plt.subplots()
ax.axis('off')  # Hide axes
legend = ax.legend(handles=legend_handles, loc='center', fontsize=fontsize, ncol=1)
plt.tight_layout()
plt.savefig('../scoring_scheme_comparison__legend.png')
plt.savefig('../scoring_scheme_comparison__legend.pdf')
plt.close()

        
