
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


from scipy.stats import pearsonr, spearmanr

from collections import Counter

import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--base_dir', type=str, default='/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/')
    parser.add_argument('--system', type=str, required=True)
    parser.add_argument('--system_name_in_csv_file', type=str, default=None)
    parser.add_argument('--target_column', type=str, required=True)
    parser.add_argument('--prediction_column', type=str, default=None)
    parser.add_argument('--model', type=str, required=True, choices=['hermes', 'proteinmpnn', 'neg_abs_diff_vdw_radius', 'blosum62', 'luksza_cross_reactivity', 'luksza_cross_reactivity_without_d'])
    parser.add_argument('--use_mt_structure', type=int, default=0, choices=[0, 1])
    parser.add_argument('--model_instance', type=str, default=None, help='Only used if model is hermes')
    parser.add_argument('--num_seq_per_target', type=int, default=None, help='Only used if model is proteinmpnn')
    parser.add_argument('--show_wt_lines', type=str, default='none', choices=['both', 'pred', 'target', 'none', 'both_from_df'])
    parser.add_argument('--wt_value_target', type=float, default=0.0)
    parser.add_argument('--wt_value_pred', type=float, default=0.0)
    args = parser.parse_args()

    system = args.system
    if args.system_name_in_csv_file is None:
        system_name_in_csv_file = system
    else:
        system_name_in_csv_file = args.system_name_in_csv_file
    target_column = args.target_column
    model = args.model
    model_instance = args.model_instance
    if model == 'hermes': assert model_instance is not None
    use_mt_structure = args.use_mt_structure
    num_seq_per_target = args.num_seq_per_target
    if model == 'proteinmpnn': assert num_seq_per_target is not None


    if model == 'hermes':
        base_dir = os.path.join(args.base_dir, f'{system}/results/{model_instance}/')
        df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-{model_instance}-use_mt_structure={use_mt_structure}.csv'))
        out_file = f'pretty_scatterplot-{system_name_in_csv_file}-{model}-{model_instance}-use_mt_structure={use_mt_structure}.png'
        title = f'{system_name_in_csv_file}\n{model} - {model_instance}\nuse_mt_structure={use_mt_structure}'
        if args.prediction_column is None:
            prediction_column = 'log_proba_mt__minus__log_proba_wt'
            ylabel = r'HERMES predictions, $\Delta logP$'
        else:
            prediction_column = args.prediction_column
            ylabel = r'HERMES predictions, ' + prediction_column
            out_file = out_file.replace('pretty_scatterplot', f'pretty_scatterplot-{prediction_column}')
        ylabel = r'HERMES predictions, $\Delta logP$'
        color = 'tab:purple'
    elif model == 'proteinmpnn':
        base_dir = os.path.join(args.base_dir, f'{system}/results/{model_instance}/zero_shot_predictions/')
        df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-num_seq_per_target={num_seq_per_target}-use_mt_structure={use_mt_structure}.csv'))
        out_file = f'pretty_scatterplot-{system_name_in_csv_file}-{model_instance}-num_seq_per_target={num_seq_per_target}-use_mt_structure={use_mt_structure}.png'
        title = f'{system_name_in_csv_file}\n{model} - num_seq_per_target={num_seq_per_target}\nuse_mt_structure={use_mt_structure}'
        if args.prediction_column is None:
            prediction_column = 'log_p_mt__minus__log_p_wt'
            ylabel = r'ProteinMPNN pred, $\Delta logP$'
        else:
            prediction_column = args.prediction_column
            ylabel = r'ProteinMPNN pred, ' + prediction_column
            out_file = out_file.replace('pretty_scatterplot', f'pretty_scatterplot-{prediction_column}')
        color = 'tab:brown'
    elif model in {'blosum62', 'neg_abs_diff_vdw_radius', 'luksza_cross_reactivity', 'luksza_cross_reactivity_without_d'}:
        print('Note: "use_mt_structure" is irrelevant with neg_abs_diff_vdw_radius substitution matrix.', file=sys.stderr)
        base_dir = os.path.join(args.base_dir, f'{system}/results/{model_instance}/')
        df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-{model_instance}-use_mt_structure=0.csv'))
        out_file = f'pretty_scatterplot-{system_name_in_csv_file}-{model_instance}.png'
        title = f'{system_name_in_csv_file}\n{model_instance}'
        prediction_column = 'substitution_matrix_score'
        ylabel = r'Substitution Matrix Score'
        color = 'tab:pink'

    # exlude nans, keep only one wildtype measurement
    is_not_nan_mask = np.logical_and(~np.isnan(df_full[target_column]), ~np.isnan(df_full[prediction_column]))
    
    if 'is_wt' in df_full.columns:
        is_wt_mask = df_full['is_wt']
        if np.sum(is_wt_mask) > 1:
            first_wt_mask = np.zeros(df_full.shape[0], dtype=bool)
            first_wt_mask[np.arange(df_full.shape[0])[is_wt_mask][0]] = True
            mask = np.logical_or(np.logical_and(is_not_nan_mask, ~is_wt_mask), first_wt_mask)
        else:
            mask = is_not_nan_mask
    else:
        print('Warning: no "is_wt" column in csv file, assuming there are no wildtype duplicates.', file=sys.stderr)
        mask = is_not_nan_mask

    predictions = df_full[prediction_column][mask]
    targets = df_full[target_column][mask]

    is_bad_value_mask = np.zeros(targets.shape[0], dtype=bool) # placeholder for potential future use

    pr_fil, pr_pval_fil = pearsonr(targets[~is_bad_value_mask], predictions[~is_bad_value_mask])
    sr_fil, sr_pval_fil = spearmanr(targets[~is_bad_value_mask], predictions[~is_bad_value_mask])
    pr, pr_pval = pearsonr(targets, predictions)
    sr, sr_pval = spearmanr(targets, predictions)
    num = targets.shape[0]

    plt.figure(figsize=(3.5, 3.5))

    if args.show_wt_lines in {'both', 'target'}:
        plt.axvline(x=args.wt_value_target, c='k', alpha=1.0, ls='--')

    if args.show_wt_lines in {'both', 'pred'}:
        plt.axhline(y=args.wt_value_pred, c='k', alpha=1.0, ls='--')
    
    if args.show_wt_lines == 'both_from_df':
        wt_values = df_full[df_full['is_wt'].astype(bool)][[target_column, prediction_column]].values
        print(wt_values)
        for wt_value_target, wt_value_pred in wt_values:
            plt.axvline(x=wt_value_target, c='k', alpha=1.0, ls='--')
            plt.axhline(y=wt_value_pred, c='k', alpha=1.0, ls='--')

    print(predictions.shape, targets.shape)

    plt.scatter(targets[is_bad_value_mask], predictions[is_bad_value_mask], c='tab:grey', alpha=0.2)
    plt.scatter(targets[~is_bad_value_mask], predictions[~is_bad_value_mask], c=color, alpha=0.5)

    plt.xlabel(args.target_column, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.title('Pearson corr. = %.2f (%.2f)\nSpearman corr. = %.2f (%.2f)' % (pr, pr_fil, sr, sr_fil), fontsize=12)
    # plt.legend(loc='upper left')
    plt.tight_layout()
    plt.savefig(os.path.join(base_dir, out_file), dpi=300)
    plt.savefig(os.path.join(base_dir, out_file.replace('.png', '.pdf')), dpi=300)
    plt.close()


