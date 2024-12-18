
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

import argparse

IND_TO_AA = 'CSTAGPDEQNHRKMILVWYF'
AA_TO_IND = {aa: i for i, aa in enumerate(IND_TO_AA)}

SYSTEM_TO_WT_TO_SHOW_COLUMNS = {
    'nyeso': None,
    'tax': None,
    'mart': None,
    'hsiue_et_al': 'is_wt',
    'mskcc': 'is_wt'
}

SYSTEM_TO_TARGET_COLUMN = {
    'nyeso': '-log10(Kd)',
    'tax': '-log10(Kd)',
    'mart': '-log10(Kd)',
    'hsiue_et_al': 'IFN_gamma (pg/ml)',
    'mskcc': '- delta log_ec50_M'
}

# SYSTEM_TO_PRETTY_TARGET = {
#     'nyeso': '$-\\text{log}_{10}(K_d)$',
#     'tax': '$-\\text{log}_{10}(K_d)$',
#     'mart': '$-\\text{log}_{10}(K_d)$',
#     'hsiue_et_al': 'IFN_gamma (pg/ml)',
#     'mskcc': '$- \Delta \\text{log(EC50)}$'
# }

SYSTEM_TO_PRETTY_TARGET = {
    'nyeso': '$-\\text{log}_{10}(K_d)$',
    'tax': '$-\\text{log}_{10}(K_d)$',
    'mart': '$-\\text{log}_{10}(K_d)$',
    'hsiue_et_al': 'IFN_gamma (pg/ml)',
    'mskcc': '$- \Delta \\text{log(EC50)}$'
}

SYSTEM_TO_MARKER_SIZE = {
    'nyeso': 60,
    'tax': 60,
    'mart': 60,
    'hsiue_et_al': 30,
    'mskcc': 30
}

SYSTEM_TO_ALPHA = {
    'nyeso': 0.6,
    'tax': 0.6,
    'mart': 0.6,
    'hsiue_et_al': 0.6,
    'mskcc': 0.6
}

MOTIF_TO_MARKER = {
    'GIG': 'd',
    'DRG': 'x'
}

SYSTEM_CSV_TO_TITLE = {
    'nyeso_peptide_kd_closest': 'NY-ESO-1',
    'tax_peptide_kd_closest': 'Tax',
    'mart_peptide_kd_closest': 'MART-1 (right motifs)',
    'mart_peptide_kd_3qdg': 'MART-1 (3qdg, GIG motif)',
    'mart_peptide_kd_6am5': 'MART-1 (6am5, GIG motif)',
    'mart_peptide_kd_6amu': 'MART-1 (6amu, DRG motif)',
    'hsiue_et_al_H2_sat_mut': 'H2-scDb on p53$^\\text{R175H}$',
    'mskcc_tcr1_ec50_sat_mut_af3': 'TCR1 on NLVPMVATV',
    'mskcc_tcr2_ec50_sat_mut_af3': 'TCR2 on NLVPMVATV (AF3)',
    'mskcc_tcr3_ec50_sat_mut_af3': 'TCR3 on NLVPMVATV (AF3)',
    'mskcc_tcr4_ec50_sat_mut_af3': 'TCR4 on IMDQVPFSV (AF3)',
    'mskcc_tcr5_ec50_sat_mut_af3': 'TCR5 on IMDQVPFSV (AF3)',
    'mskcc_tcr6_ec50_sat_mut_af3': 'TCR6 on IMDQVPFSV (AF3)',
    'mskcc_tcr7_ec50_sat_mut_af3': 'TCR7 on GRLKALCQR (AF3)',
}

'''
Scatterplot for chosen HERMES model
'''



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_version', type=str, required=True)
    parser.add_argument('--with_relaxation', type=int, choices=[0, 1], required=True)
    parser.add_argument('--system', type=str, required=True)
    parser.add_argument('--system_name_in_csv_file', type=str, required=True)

    args = parser.parse_args()

    target_column = SYSTEM_TO_TARGET_COLUMN[args.system]
    prediction_column = 'pnE'
    title = SYSTEM_CSV_TO_TITLE[args.system_name_in_csv_file]
    xlabel = SYSTEM_TO_PRETTY_TARGET[args.system]
    ylabel = r'$E(\sigma_p ; \text{TCR}, \text{MHC})$'
    color = 'tab:purple'
    s = SYSTEM_TO_MARKER_SIZE[args.system]
    alpha = SYSTEM_TO_ALPHA[args.system]
    
    if args.model_version == 'hermes_py_000' and args.with_relaxation:
        model_pretty_name = 'HERMES-relaxed 0.00'
    elif args.model_version == 'hermes_py_000' and not args.with_relaxation:
        model_pretty_name = 'HERMES-fixed 0.00'
    elif args.model_version == 'hermes_py_050' and args.with_relaxation:
        model_pretty_name = 'HERMES-relaxed 0.50'
    elif args.model_version == 'hermes_py_050' and not args.with_relaxation:
        model_pretty_name = 'HERMES-fixed 0.50'
    else:
        raise ValueError(f'Unknown model version: {args.model_version}')

    with_relaxation_str = '_with_relaxation' if args.with_relaxation else ''

    df = pd.read_csv(f'../mutation_effects/{args.system}/results/{args.model_version}/{args.system_name_in_csv_file}{with_relaxation_str}-{args.model_version}-use_mt_structure=0.csv')

    # exlude nans, keep only one wildtype measurement
    is_not_nan_mask = np.logical_and(~np.isnan(df[target_column]), ~np.isnan(df[prediction_column]))
    
    if 'is_wt' in df.columns:
        is_wt_mask = df['is_wt'].astype(bool)
        first_wt_mask = np.zeros(df.shape[0], dtype=bool)
        first_wt_mask[np.arange(df.shape[0])[is_wt_mask][0]] = True
        mask = np.logical_or(np.logical_and(is_not_nan_mask, ~is_wt_mask), first_wt_mask)
    else:
        print('Warning: no "is_wt" column in csv file, assuming there are no wildtype duplicates.', file=sys.stderr)
        mask = is_not_nan_mask

    df = df[mask]

    predictions = df[prediction_column]
    targets = df[target_column]

    if 'motif' in df.columns:
        markers = np.array([MOTIF_TO_MARKER[motif] for motif in df['motif'].values])
    else:
        markers = np.array(['o'] * df.shape[0])

    pr, pr_pval = pearsonr(targets, predictions)
    sr, sr_pval = spearmanr(targets, predictions)
    num = targets.shape[0]

    plt.figure(figsize=(3.5, 3.5))

    
    if SYSTEM_TO_WT_TO_SHOW_COLUMNS[args.system] is not None:
        wt_values = df[df[SYSTEM_TO_WT_TO_SHOW_COLUMNS[args.system]].astype(bool)][[target_column, prediction_column]].values
        for wt_value_target, wt_value_pred in wt_values:
            plt.axvline(x=wt_value_target, c='k', alpha=1.0, ls='--', lw=0.85)
            plt.axhline(y=wt_value_pred, c='k', alpha=1.0, ls='--', lw=0.85)

    unique_markers = np.unique(markers)
    for marker in unique_markers:
        mask = markers == marker
        plt.scatter(targets[mask], predictions[mask], c=color, marker=marker, alpha=alpha, s=s)
    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)

    plt.title(title, fontsize=15)
    # plt.title(title+'\n'+f'Spearman r = {sr:.2f} (pv {sr_pval:.1e})', fontsize=12)

    # put text of correlation coefficient
    plt.text(0.03, 0.03, f'Spearman r = {sr:.2f}\np-val = {sr_pval:.1e}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='left')

    # plt.legend(loc='upper left')

    plt.tight_layout()
    outfile = f'../mutation_effects/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-scatterplot.png'
    plt.savefig(outfile, dpi=300)
    plt.savefig(outfile.replace('.png', '.pdf'), dpi=300)
    plt.close()

    if len(unique_markers) > 1:
        # make legend in standlone figure
        plt.figure(figsize=(1.5, 1.5))
        for motif in MOTIF_TO_MARKER:
            marker = MOTIF_TO_MARKER[motif]
            plt.scatter([], [], c=color, marker=marker, alpha=alpha, s=s, label=motif)
        plt.legend(loc='center')
        plt.axis('off')
        plt.tight_layout()
        outfile = f'../mutation_effects/{args.system}/plots/motif_legend.png'
        plt.savefig(outfile, dpi=300)
        plt.savefig(outfile.replace('.png', '.pdf'), dpi=300)
        plt.close()
    
    # now also make two heatmaps! one with the predictions, one with the targets
    plt.figure(figsize=(3.5, 3.5))

    sequences = df['sequence']
    predictions = df[prediction_column]
    targets = df[target_column]

    heatmap = np.full((len(sequences[0]), 20), np.nan)
    

