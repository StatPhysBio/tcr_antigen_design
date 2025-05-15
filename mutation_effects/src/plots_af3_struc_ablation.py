
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from scipy.stats import spearmanr
import argparse

from mutation_effects_comparison_plots import make_pretty_name, get_long_prediction_column_name, MUT_EFFECTS_DIR, SYSTEM_TO_TARGET_COLUMN

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

MODEL_TO_COLOR = {
    ('hermes_py_000', 'pE-relaxed'): orange,
    ('hermes_py_050', 'pE-relaxed'): orange_light,
    ('hermes_py_000', 'pE-fixed'): red,
    ('hermes_py_050', 'pE-fixed'): red_light,
    ('proteinmpnn_v_48_002_full_pep_masked', 'pnlogp-fixed'): blue,
    ('proteinmpnn_v_48_020_full_pep_masked', 'pnlogp-fixed'): blue_light,
}

SYSTEM_TO_PRETTY_NAME = {
    ('nyeso', None): '1G4 TCR',
    ('tax', None): 'A6 TCR',
    ('hsiue_et_al', None): 'H2-scDb (Fab)',
    ('mskcc', 1): 'TCR1',
    ('mskcc', 2): 'TCR2',
    ('mskcc', 3): 'TCR3',
    ('mskcc', 4): 'TCR4',
    ('mskcc', 5): 'TCR5',
    ('mskcc', 6): 'TCR6',
    ('mskcc', 7): 'TCR7',
}

STRUC_TO_PRETTY_NAME = {
    'true struc': 'Experim.',
    'af3 yes template': 'AF3\nw/ templ.',
    'af3 no template': 'AF3\nw/out templ.',
}

def system_to_csv_names(system, tcr=None):

    if system == 'nyeso':
        return {'true struc': 'nyeso_peptide_kd_closest', 'af3 yes template': 'nyeso_peptide_kd_closest_af3_yes_template', 'af3 no template': 'nyeso_peptide_kd_closest_af3_no_template'}
    elif system == 'tax':
        return {'true struc': 'tax_peptide_kd_closest', 'af3 yes template': 'tax_peptide_kd_closest_af3_yes_template', 'af3 no template': 'tax_peptide_kd_closest_af3_no_template'}
    elif system == 'hsiue_et_al':
        return {'true struc': 'hsiue_et_al_H2_sat_mut', 'af3 yes template': 'hsiue_et_al_H2_sat_mut_af3_yes_template', 'af3 no template': 'hsiue_et_al_H2_sat_mut_af3_no_template'}
    elif system == 'mskcc':
        if tcr == 1:
            return {'true struc': 'mskcc_tcr1_ec50_sat_mut_af3', 'af3 yes template': 'mskcc_tcr1_ec50_sat_mut_af3_yes_template', 'af3 no template': 'mskcc_tcr1_ec50_sat_mut_af3_yes_template'}
        elif tcr >= 2 and tcr <= 7:
            return {'true struc': None, 'af3 yes template': f'mskcc_tcr{tcr}_ec50_sat_mut_af3', 'af3 no template': f'mskcc_tcr{tcr}_ec50_sat_mut_af3_no_template'}
        else:
            raise ValueError(f'Invalid TCR value: {tcr}')
    else:
        raise ValueError(f'Invalid system: {system}')


def get_model_specific_parameters(model_instance, prediction_column_short, system, tcr=None):

    use_mt_structure = 0 # always 0, there for legacy reasons
    num_seq_per_target = 10 # always 10, there for legacy reasons

    color = MODEL_TO_COLOR[(model_instance, prediction_column_short)]

    struc_to_df = {}

    if 'hermes' in model_instance:
        base_dir = f'{MUT_EFFECTS_DIR}/{system}/results/{model_instance}/'

        struc_to_csvname = system_to_csv_names(system, tcr)
        for struc in struc_to_csvname.keys():
            if struc_to_csvname[struc] is None:
                struc_to_df[struc] = None
                continue
            prediction_column, system_name_in_csv_file = get_long_prediction_column_name(model_instance, prediction_column_short, struc_to_csvname[struc])
            try:
                struc_to_df[struc] = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-{model_instance}-use_mt_structure={use_mt_structure}.csv'))
            except FileNotFoundError:
                struc_to_df[struc] = None

    elif 'proteinmpnn' in model_instance:
        base_dir = f'{MUT_EFFECTS_DIR}/{system}/results/{model_instance}/zero_shot_predictions/'

        struc_to_csvname = system_to_csv_names(system, tcr)
        for struc in struc_to_csvname.keys():
            if struc_to_csvname[struc] is None:
                struc_to_df[struc] = None
                continue
            prediction_column, system_name_in_csv_file = get_long_prediction_column_name(model_instance, prediction_column_short, struc_to_csvname[struc])
            try:
                struc_to_df[struc] = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-num_seq_per_target={num_seq_per_target}-use_mt_structure={use_mt_structure}.csv'))
            except FileNotFoundError:
                struc_to_df[struc] = None

    else:

        raise ValueError(f'Invalid model_instance: {model_instance}')
    
    return struc_to_df, prediction_column, color


def get_spearmanr(df, prediction_column, target_column):

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

    
    if 'is_reliable' in df.columns:
        is_reliable_mask = np.logical_and(df['is_reliable'].values, mask)
        # is_not_reliable_mask = np.logical_and(~df['is_reliable'].values, mask)

        predictions_rel = df[prediction_column][is_reliable_mask]
        targets_rel = df[target_column][is_reliable_mask]

        # pr, pr_pval = pearsonr(targets_rel, predictions_rel)
        sr, sr_pval = spearmanr(targets_rel, predictions_rel)

    else:
        predictions = df[prediction_column][mask]
        targets = df[target_column][mask]

        # pr, pr_pval = pearsonr(targets, predictions)
        sr, sr_pval = spearmanr(targets, predictions)
    
    return sr, sr_pval



if __name__ == '__main__':

    system_and_tcr = [('nyeso', None), ('tax', None), ('hsiue_et_al', None), ('mskcc', 1), ('mskcc', 2), ('mskcc', 3), ('mskcc', 4), ('mskcc', 5), ('mskcc', 6), ('mskcc', 7)]

    models = ['hermes_py_000', 'hermes_py_050', 'hermes_py_000', 'hermes_py_050', 'proteinmpnn_v_48_002_full_pep_masked', 'proteinmpnn_v_48_020_full_pep_masked']
    pred_columns_short = ['pE-fixed', 'pE-fixed', 'pE-relaxed', 'pE-relaxed', 'pnlogp-fixed', 'pnlogp-fixed']

    fontsize = 15

    ncols = 2
    nrows = 5
    colsize = 3
    rowsize = 2.5
    fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(colsize*ncols, rowsize*nrows), dpi=300, sharex=True, sharey=True)

    ylim = (np.inf, -np.inf)

    for ax_i, (system, tcr) in enumerate(system_and_tcr):
        ax = axs[ax_i//ncols, ax_i%ncols]
        
        title = SYSTEM_TO_PRETTY_NAME[(system, tcr)]

        for j, (model_instance, pred_column_short) in enumerate(zip(models, pred_columns_short)):

            struc_to_df, prediction_column, color = get_model_specific_parameters(model_instance, pred_column_short, system, tcr)

            xticks = []
            srs = []
            srs_pvals = []
            for struc_i, struc in enumerate(struc_to_df.keys()):
                df = struc_to_df[struc]

                if df is None:
                    sr, sr_pval = np.nan, np.nan
                else:
                    sr, sr_pval = get_spearmanr(df, prediction_column, SYSTEM_TO_TARGET_COLUMN[system])

                xticks.append(STRUC_TO_PRETTY_NAME[struc])
                srs.append(sr)
                srs_pvals.append(sr_pval)
            
            ax.plot(np.arange(len(xticks)), srs, marker='o', color=color, alpha=0.8, linewidth=2.5)

            ax.set_xticks(np.arange(len(xticks)))
            ax.set_xticklabels(xticks, fontsize=fontsize)

            ax.tick_params(axis='x', labelrotation=60, labelsize=fontsize)
            
            if ax_i%ncols == 0:
                ax.set_ylabel('Spearman r', fontsize=fontsize)
            
            ax.tick_params(axis='y', labelsize=fontsize-2)

            ax.set_title(title, fontsize=fontsize+2)
            ax.axhline(0, color='black', linestyle='--', linewidth=1.0)

            ax.grid(color='gray', linestyle='--', alpha=0.5)

            ylim = (min(ylim[0], np.nanmin(srs)-0.075), max(ylim[1], np.nanmax(srs)+0.075))
    
    # samw ylim but qwithout ahrey=True because I want the yticklabels to be displayed
    for ax in axs.flatten():
        ax.set_ylim(ylim)

    plt.tight_layout()
    plt.savefig(os.path.join(MUT_EFFECTS_DIR, 'plots_af3_struc_ablation.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(MUT_EFFECTS_DIR, 'plots_af3_struc_ablation.png'), bbox_inches='tight')
    plt.close()


    MODEL_TO_COLOR = {
        'hermes_fixed_000': red,
        'hermes_fixed_050': red_light,
        'hermes_relaxed_000': orange,
        'hermes_relaxed_050': orange_light,
        'proteinmpnn_002': blue,
        'proteinmpnn_020': blue_light,
    }

    MODEL_TO_PRETTY_NAME = {
        'hermes_fixed_000': 'HERMES-$fixed$ 0.00',
        'hermes_fixed_050': 'HERMES-$fixed$ 0.50',
        'hermes_relaxed_000': 'HERMES-$relaxed$ 0.00',
        'hermes_relaxed_050': 'HERMES-$relaxed$ 0.50',
        'proteinmpnn_002': 'ProteinMPNN 0.02 Pep',
        'proteinmpnn_020': 'ProteinMPNN 0.20 Pep',
    }

    ## make legend
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    handles = []
    for model in MODEL_TO_COLOR.keys():
        color = MODEL_TO_COLOR[model]
        model_pretty_name = MODEL_TO_PRETTY_NAME[model]
        handles.append(Patch(facecolor=color, edgecolor=color, alpha=0.8, label=model_pretty_name))
    plt.figure(figsize=(5, 5))
    plt.legend(handles=handles, loc='center', fontsize=14)
    plt.axis('off')
    plt.savefig(f'../plots_af3_struc_ablation__legend.png')
    plt.savefig(f'../plots_af3_struc_ablation__legend.pdf')
    plt.close()

