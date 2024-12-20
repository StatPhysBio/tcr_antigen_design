
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr
import json

import logging
logging.basicConfig(level=logging.ERROR)

import argparse


## not needed for the other systems
SYSTEM_TO_PDBS = {
    'nyeso': ['2bnr', '2bnq'],
    'tax': ['1ao7', '1qse', '1qsf'],
    'mart': ['3qdg', '6am5', '6amu']
}

# None means that the wildtype won't be shown
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
    'mskcc': '-log10(EC50)'
}

def get_long_prediction_column_name(model_instance, prediction_column_short, system_name_in_csv_file):
    if 'hermes' in model_instance:
        if prediction_column_short == 'delta_log_p':
            return 'log_proba_mt__minus__log_proba_wt', system_name_in_csv_file
        elif prediction_column_short == 'pE-relaxed':
            return 'pnE', system_name_in_csv_file + '_with_relaxation'
        elif prediction_column_short == 'pE-fixed':
            return 'pnE', system_name_in_csv_file
        else:
            raise ValueError(f'Unknown prediction_column_short: {prediction_column_short}')
        
    elif 'proteinmpnn' in model_instance:
        if prediction_column_short == 'delta_log_p':
            return 'log_p_mt__minus__log_p_wt', system_name_in_csv_file
        elif prediction_column_short == 'pnlogp-fixed':
            return 'pnlogp', system_name_in_csv_file
        else:
            raise ValueError(f'Unknown prediction_column_short: {prediction_column_short}')
        
    elif 'tcrdock' in model_instance:
        system_name_in_csv_file = system_name_in_csv_file.replace('_averaged', '')
        system_name_in_csv_file = system_name_in_csv_file.replace('_closest', '')
        if prediction_column_short == 'neg_pae':
            return 'neg pmhc_tcr_pae', system_name_in_csv_file
        else:
            raise ValueError(f'Unknown prediction_column_short: {prediction_column_short}')
    
    elif 'blosum62' in model_instance:
        if prediction_column_short == 'sub_score':
            return 'substitution_matrix_score', system_name_in_csv_file
        else:
            raise ValueError(f'Unknown prediction_column_short: {prediction_column_short}')
    
    else:
        raise ValueError(f'Unknown model_instance: {model_instance}')


SHORT_PREDICTION_NAME_TO_PRETTY_NAME = {
    'pE-fixed': '$pE$ on fixed structure',
    'pE-relaxed': 'avg. $pE$ on mutated+relaxed structures',
    'neg_pae': '-PAE',
    'sub_score': 'Substitution Matrix Score',
    'pnlogp-fixed': 'log P on fixed structure',
    'delta_log_p': '$pE$ on fixed structure'
}

def short_prediction_name_to_pretty_name(short_prediction_name):
    if short_prediction_name == 'pE-fixed':
        return '$pE$ on fixed structure'
    elif short_prediction_name == 'pE-relaxed':
        return 'avg. $pE$ on\nmutated+relaxed structures'
    elif short_prediction_name == 'neg_pae':
        return '-PAE'
    elif short_prediction_name == 'sub_score':
        return 'Substitution Matrix Score'
    elif short_prediction_name == 'pnlogp-fixed':
        return 'log P on fixed structure'
    elif short_prediction_name == 'delta_log_p':
        return '$pE$ on fixed structure'

MODEL_INSTANCE_TO_PRETTY_NAME = {
    'hermes_py_000': 'HERMES 0.00',
    'hermes_py_050': 'HERMES 0.50',
    'hermes_py_000_ft_skempi_no_tcrpmhc_ddg_bi': 'HERMES 0.00 + Skempi FT',
    'hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi': 'HERMES 0.50 + Skempi FT',
    'proteinmpnn_v_48_002': 'ProteinMPNN 0.02',
    'proteinmpnn_v_48_030': 'ProteinMPNN 0.30',
    'tcrdock': 'TCRdock',
    'tcrdock_no_nearby_templates': 'TCRdock benchmark',
    'blosum62': 'BLOSUM62'
}

def make_pretty_name(model, pred_col):
    if model == 'hermes_py_000' and pred_col in {'delta_log_p', 'pE-fixed'}:
        return 'HERMES-fixed 0.00'
    elif model == 'hermes_py_050' and pred_col in {'delta_log_p', 'pE-fixed'}:
        return 'HERMES-fixed 0.50'
    elif model == 'hermes_py_000' and pred_col == 'pE-relaxed':
        return 'HERMES-relaxed 0.00'
    elif model == 'hermes_py_050' and pred_col == 'pE-relaxed':
        return 'HERMES-relaxed 0.50'
    elif model == 'proteinmpnn_v_48_002' and pred_col in {'delta_log_p', 'pnlogp-fixed'}:
        return 'ProteinMPNN 0.02'
    elif model == 'proteinmpnn_v_48_030' and pred_col in {'delta_log_p', 'pnlogp-fixed'}:
        return 'ProteinMPNN 0.30'
    elif model == 'tcrdock' and pred_col == 'neg_pae':
        return 'TCRdock'
    elif model == 'tcrdock_no_nearby_templates' and pred_col == 'neg_pae':
        return 'TCRdock benchmark'
    elif model == 'blosum62' and pred_col == 'sub_score':
        return 'BLOSUM62'
    else:
        raise ValueError(f'Unknown model and pred_col: {model}, {pred_col}')


def get_model_specific_parameters(model_instance, prediction_column_short, args):

    system = args.system
    system_name_in_csv_file = args.system_name_in_csv_file
    use_mt_structure = 0 # always 0, there for legagy reasons
    num_seq_per_target = 10 # always 10, there for legagy reasons

    prediction_column, system_name_in_csv_file = get_long_prediction_column_name(model_instance, prediction_column_short, system_name_in_csv_file)
    title = make_pretty_name(model_instance, prediction_column_short)

    if 'hermes' in model_instance:
        base_dir = f'../mutation_effects/{system}/results/{model_instance}/'
        color = 'tab:purple'

        df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-{model_instance}-use_mt_structure={use_mt_structure}.csv'))

    elif 'proteinmpnn' in model_instance:
        base_dir = f'../mutation_effects/{system}/results/{model_instance}/zero_shot_predictions/'

        color = 'tab:brown'

        df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-num_seq_per_target={num_seq_per_target}-use_mt_structure={use_mt_structure}.csv'))

    elif model_instance == 'tcrdock':
        # print('Note: "use_mt_structure" is irrelevant with tcrdock model.', file=sys.stderr)
        base_dir = f'../mutation_effects/{system}/results/{model_instance}/'
        try:
            df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}_w_pae_filtered.tsv'), sep='\t')
        except:
            df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}_w_pae.tsv'), sep='\t')
        df_full['neg pmhc_tcr_pae'] = -df_full['pmhc_tcr_pae']
        color = 'tab:cyan'
    
    elif model_instance == 'tcrdock_no_nearby_templates':
        # print('Note: "use_mt_structure" is irrelevant with tcrdock model.', file=sys.stderr)
        base_dir = f'../mutation_effects/{system}/results/{model_instance}/'
        try:
            df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}_no_nearby_templates_w_pae_filtered.tsv'), sep='\t')
        except:
            df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}_no_nearby_templates_w_pae.tsv'), sep='\t')
        df_full['neg pmhc_tcr_pae'] = -df_full['pmhc_tcr_pae']
        color = 'tab:cyan'

    # elif model == 'neg_abs_diff_vdw_radius':
    #     print('Note: "use_mt_structure" is irrelevant with neg_abs_diff_vdw_radius substitution matrix.', file=sys.stderr)
    #     base_dir = os.path.join(args.base_dir, f'{system}/results/{model_instance}/')
    #     df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-{model_instance}-use_mt_structure=0.csv'))
    #     out_file = f'pretty_scatterplot-{system_name_in_csv_file}-{model_instance}.png'
    #     title = f'{system_name_in_csv_file}\n{model_instance}'
    #     prediction_column = 'substitution_matrix_score'
    #     ylabel = r'Substitution Matrix Score'
    #     color = 'tab:red'

    elif 'blosum62' in model_instance:
        # print('Note: "use_mt_structure" is irrelevant with neg_abs_diff_vdw_radius substitution matrix.', file=sys.stderr)
        base_dir = f'../mutation_effects/{system}/results/{model_instance}/'

        if 'average' in system_name_in_csv_file:
            df_list = [pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file.replace("averaged", pdb)}-{model_instance}-use_mt_structure={use_mt_structure}.csv')) for pdb in SYSTEM_TO_PDBS[system]]
            # df_full is same as df_list but prediction column is the average of all predictions
            df_full = df_list[0].copy()
            for df in df_list[1:]:
                df_full[prediction_column] += df[prediction_column]
            df_full[prediction_column] /= len(df_list)
        else:
            df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-{model_instance}-use_mt_structure=0.csv'))

        color = 'tab:pink'
    
    else:
        raise ValueError(f'Unknown model_instance: {model_instance}')
    
    return df_full, title, prediction_column, color


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--system', type=str, required=True)
    parser.add_argument('--system_name_in_csv_file', type=str, required=True)
    args = parser.parse_args()

    # get pretty target column
    if args.system == 'nyeso':
        target_column_pretty = r'$-log_{10}(K_d)$'
    elif args.system == 'tax':
        target_column_pretty = r'$-log_{10}(K_d)$'
    elif args.system == 'mart':
        target_column_pretty = r'$-log_{10}(K_d)$'
    elif args.system == 'hsiue_et_al':
        target_column_pretty = 'IFN_gamma (pg/ml)'
    elif args.system == 'mskcc':
        target_column_pretty = r'$-\Delta log_{EC50}$'
    else:
        raise ValueError(f'Unknown system: {args.system}')

    # make the list of models and the figure shape based on the system
    if 'hsiue' in args.system:
        models = ['hermes_py_000', 'hermes_py_050', 'hermes_py_000', 'hermes_py_050', 'proteinmpnn_v_48_002', 'proteinmpnn_v_48_030', 'blosum62']
        prediction_columns = ['pE-fixed', 'pE-fixed', 'pE-relaxed', 'pE-relaxed', 'pnlogp-fixed', 'pnlogp-fixed', 'sub_score']
        num_rows = 4
        num_cols = 2
    elif 'mskcc' in args.system:
        models = ['hermes_py_000', 'hermes_py_050', 'hermes_py_000', 'hermes_py_050', 'proteinmpnn_v_48_002', 'proteinmpnn_v_48_030', 'tcrdock', 'tcrdock_no_nearby_templates', 'blosum62']
        prediction_columns = ['pE-fixed', 'pE-fixed', 'pE-relaxed', 'pE-relaxed', 'pnlogp-fixed', 'pnlogp-fixed', 'neg_pae', 'neg_pae', 'sub_score']
        num_rows = 5
        num_cols = 2
    else:
        models = ['hermes_py_000', 'hermes_py_050', 'hermes_py_000', 'hermes_py_050', 'proteinmpnn_v_48_002', 'proteinmpnn_v_48_030', 'tcrdock', 'tcrdock_no_nearby_templates', 'blosum62']
        prediction_columns = ['pE-fixed', 'pE-fixed', 'pE-relaxed', 'pE-relaxed', 'pnlogp-fixed', 'pnlogp-fixed', 'neg_pae', 'neg_pae', 'sub_score']
        num_rows = 5
        num_cols = 2
    
    target_column = SYSTEM_TO_TARGET_COLUMN[args.system]
    # target_column_pretty = system_to_target_column_pretty(args.system)
    

    ## scatterplots
    
    rowsize = 4
    colsize = 4
    fontsize_base = 18
    fig, axs = plt.subplots(num_rows, num_cols, figsize=(num_cols*colsize, num_rows*rowsize))

    metrics = {'Pr': [], 'Sr': [], 'Pr_pval': [], 'Sr_pval': [], 'AUROC reliable': [], 'AUROC all': []}
    colors = []
    model_names_pretty = []

    for i, (model_instance, prediction_column_short) in enumerate(zip(models, prediction_columns)):
        row = i // num_cols
        col = i % num_cols
        ax = axs[row, col]

        try:
            df_full, title, prediction_column, color = get_model_specific_parameters(model_instance, prediction_column_short, args)
        except Exception as e:
            print(f'Warning: Could not find {model_instance} for {args.system} system and prediction value {prediction_column_short}.', file=sys.stderr)
            # raise e # uncomment when debugging
            continue

        # get pretty prediction column
        prediction_column_pretty = SHORT_PREDICTION_NAME_TO_PRETTY_NAME[prediction_column_short]

        # exlude nans, keep only one wildtype measurement
        is_not_nan_mask = np.logical_and(~np.isnan(df_full[target_column]), ~np.isnan(df_full[prediction_column]))
        if 'is_wt' in df_full.columns:
            is_wt_mask = df_full['is_wt'].astype(bool)
            first_wt_mask = np.zeros(df_full.shape[0], dtype=bool)
            first_wt_mask[np.arange(df_full.shape[0])[is_wt_mask][0]] = True
            mask = np.logical_or(np.logical_and(is_not_nan_mask, ~is_wt_mask), first_wt_mask)
        else:
            print('Warning: no "is_wt" column in csv file, assuming there are no wildtype duplicates.', file=sys.stderr)
            mask = is_not_nan_mask

        if SYSTEM_TO_WT_TO_SHOW_COLUMNS[args.system] is not None:
            wt_values = df_full[df_full[SYSTEM_TO_WT_TO_SHOW_COLUMNS[args.system]].astype(bool)][[target_column, prediction_column]].values
            for wt_value_target, wt_value_pred in wt_values:
                ax.axvline(x=wt_value_target, c='k', alpha=1.0, ls='--')
                ax.axhline(y=wt_value_pred, c='k', alpha=1.0, ls='--')
        
        # if system is mskcc (luksza data), consider only reliable measurements for correlations, and also compute AUROC with and without gray points (threshold is value above vs below that of wildtype)
        if args.system == 'mskcc':

            wt_target_value = wt_values[0][0] # there's only one anyways

            predictions_all = df_full[prediction_column].values
            targets_all = df_full[target_column].values

            is_reliable_mask = np.logical_and(df_full['is_reliable'].values, mask)
            is_not_reliable_mask = np.logical_and(~df_full['is_reliable'].values, mask)

            predictions_rel = df_full[prediction_column][is_reliable_mask]
            targets_rel = df_full[target_column][is_reliable_mask]

            prediction_nonrel = df_full[prediction_column][is_not_reliable_mask]
            target_nonrel = df_full[target_column][is_not_reliable_mask]

            pr, pr_pval = pearsonr(targets_rel, predictions_rel)
            sr, sr_pval = spearmanr(targets_rel, predictions_rel)

            pr_all, pr_pval_all = pearsonr(targets_all, predictions_all)
            sr_all, sr_pval_all = spearmanr(targets_all, predictions_all)

            from sklearn.metrics import roc_auc_score
            auroc_rel = roc_auc_score(targets_rel > wt_target_value, predictions_rel)
            auroc_all = roc_auc_score(targets_all > wt_target_value, predictions_all)

            metrics['Pr'].append(pr)
            metrics['Sr'].append(sr)
            metrics['Pr_pval'].append(pr_pval)
            metrics['Sr_pval'].append(sr_pval)
            metrics['AUROC reliable'].append(auroc_rel)
            metrics['AUROC all'].append(auroc_all)

            colors.append(color)
            model_names_pretty.append(make_pretty_name(model_instance, prediction_column_short))

            ax.scatter(targets_rel, predictions_rel, color=color, s=80, alpha=0.5)
            ax.scatter(target_nonrel, prediction_nonrel, color='tab:gray', s=80, alpha=0.5)

        else:
            predictions = df_full[prediction_column][mask]
            targets = df_full[target_column][mask]

            pr, pr_pval = pearsonr(targets, predictions)
            sr, sr_pval = spearmanr(targets, predictions)

            metrics['Pr'].append(pr)
            metrics['Sr'].append(sr)
            metrics['Pr_pval'].append(pr_pval)
            metrics['Sr_pval'].append(sr_pval)
            colors.append(color)
            model_names_pretty.append(make_pretty_name(model_instance, prediction_column_short))

            ax.scatter(targets, predictions, color=color, s=80, alpha=0.5)

        ax.set_title(title, fontsize=fontsize_base+2)
        ax.set_xlabel(target_column_pretty, fontsize=fontsize_base)
        ax.set_ylabel(prediction_column_pretty, fontsize=fontsize_base)
        ax.text(0.05, 0.95, f'Pr={pr:.2f}, p={pr_pval:.2f}\nSr={sr:.2f}, p={sr_pval:.2f}', transform=ax.transAxes, fontsize=fontsize_base, verticalalignment='top')

        ax.tick_params(axis='both', which='major', labelsize=fontsize_base-2)

    os.makedirs(f'../mutation_effects/{args.system}/plots', exist_ok=True)
    plt.tight_layout()
    plt.savefig(f'../mutation_effects/{args.system}/plots/{args.system_name_in_csv_file}_scatterplots.png')
    plt.savefig(f'../mutation_effects/{args.system}/plots/{args.system_name_in_csv_file}_scatterplots.pdf')
    plt.close()
    

    ## barplots, horizontal so names are easier to read

    # invert everything
    for metric in metrics:
        metrics[metric] = metrics[metric][::-1]
    colors = colors[::-1]
    model_names_pretty = model_names_pretty[::-1]

    x = np.arange(len(model_names_pretty))
    width = 0.50

    data = {
        'Pearson r': {},
        'Spearman r': {},
        'AUROC reliable': {},
        'AUROC all': {}
    }

    fig, axs = plt.subplots(figsize=(14, num_rows*0.8), nrows=1, ncols=2, sharey=True)

    for ax, metric, metric_pretty_name in zip(axs, ['Pr', 'Sr'], ['Pearson r', 'Spearman r']):

        # barplot with the different colors
        ax.barh(x, metrics[metric], width, color=colors)

        # put values on the plot
        for i, (r, r_pval, color) in enumerate(zip(metrics[metric], metrics[f'{metric}_pval'], colors)):
            if r < 0:
                ha='right'
                offset = -0.01
            else:
                ha='left'
                offset = 0.01
            
            ax.text(r+offset, i, f'{r:.2f}', va='center', fontsize=fontsize_base-2, color='black' if r_pval <= 0.05 else 'red', ha=ha)

            data[metric_pretty_name][model_names_pretty[i]] = (r, r_pval)

        ax.set_yticks(x, model_names_pretty)

        ax.set_xlabel(metric_pretty_name, fontsize=fontsize_base)

        # change xlim to accomodate for the text
        text_offset = 0.15
        ax.set_xlim(ax.get_xlim()[0] - text_offset, ax.get_xlim()[1] + text_offset)

        ax.tick_params(axis='both', which='major', labelsize=fontsize_base-2)


    os.makedirs(f'../mutation_effects/{args.system}/plots', exist_ok=True)
    plt.tight_layout()
    plt.savefig(f'../mutation_effects/{args.system}/plots/{args.system_name_in_csv_file}_barplots.png')
    plt.savefig(f'../mutation_effects/{args.system}/plots/{args.system_name_in_csv_file}_barplots.pdf')
    plt.close()

    if args.system == 'mskcc':
        for metric in ['AUROC reliable', 'AUROC all']:
            for i in range(len(model_names_pretty)):
                data[metric][model_names_pretty[i]] = metrics[metric][i]

    # save the data
    with open(f'../mutation_effects/{args.system}/results/{args.system_name_in_csv_file}_metrics.json', 'w') as f:
        json.dump(data, f, indent=4)

    









