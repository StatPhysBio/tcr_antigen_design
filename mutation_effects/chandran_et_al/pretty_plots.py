
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
    parser.add_argument('--model', type=str, required=True, choices=['hermes', 'proteinmpnn', 'tcrdock', 'tcrdock_no_nearby_templates', 'neg_abs_diff_vdw_radius', 'blosum62'])
    parser.add_argument('--use_mt_structure', type=int, default=0, choices=[0, 1])
    parser.add_argument('--model_instance', type=str, default=None, help='Only used if model is hermes')
    parser.add_argument('--num_seq_per_target', type=int, default=None, help='Only used if model is proteinmpnn')
    parser.add_argument('--show_wt_lines', type=str, default='both', choices=['both', 'pred', 'target', 'none'])
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
        prediction_column = 'log_proba_mt__minus__log_proba_wt'
        ylabel = r'HERMES predictions, $\Delta logP$'
        color = 'tab:purple'
    elif model == 'proteinmpnn':
        base_dir = os.path.join(args.base_dir, f'{system}/results/{model_instance}/zero_shot_predictions/')
        df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-num_seq_per_target={num_seq_per_target}-use_mt_structure={use_mt_structure}.csv'))
        out_file = f'pretty_scatterplot-{system_name_in_csv_file}-{model_instance}-num_seq_per_target={num_seq_per_target}-use_mt_structure={use_mt_structure}.png'
        title = f'{system_name_in_csv_file}\n{model} - num_seq_per_target={num_seq_per_target}\nuse_mt_structure={use_mt_structure}'
        prediction_column = 'log_p_mt__minus__log_p_wt'
        ylabel = r'ProteinMPNN pred, $\Delta logP$'
        color = 'tab:brown'
    elif model == 'tcrdock':
        print('Note: "use_mt_structure" is irrelevant with tcrdock model.', file=sys.stderr)
        base_dir = os.path.join(args.base_dir, f'{system}/{model}/output/')
        df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}_w_pae.tsv'), sep='\t')
        out_file = f'pretty_scatterplot-{system_name_in_csv_file}-{model}.png'
        title = f'{system_name_in_csv_file}\n{model} - {model_instance}\nuse_mt_structure={use_mt_structure}'
        df_full['neg pmhc_tcr_pae'] = -df_full['pmhc_tcr_pae']
        prediction_column = 'neg pmhc_tcr_pae'
        ylabel = r'TCRdock predictions, neg PAE'
        color = 'tab:cyan'
    elif model == 'tcrdock_no_nearby_templates':
        print('Note: "use_mt_structure" is irrelevant with tcrdock model.', file=sys.stderr)
        base_dir = os.path.join(args.base_dir, f'{system}/tcrdock/output/')
        df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}_no_nearby_templates_w_pae.tsv'), sep='\t')
        out_file = f'pretty_scatterplot-{system_name_in_csv_file}-{model}.png'
        title = f'{system_name_in_csv_file}\n{model} - {model_instance}\nuse_mt_structure={use_mt_structure}'
        df_full['neg pmhc_tcr_pae'] = -df_full['pmhc_tcr_pae']
        prediction_column = 'neg pmhc_tcr_pae'
        ylabel = r'TCRdock predictions, neg PAE'
        color = 'tab:cyan'
    elif model == 'neg_abs_diff_vdw_radius':
        print('Note: "use_mt_structure" is irrelevant with neg_abs_diff_vdw_radius substitution matrix.', file=sys.stderr)
        base_dir = os.path.join(args.base_dir, f'{system}/results/{model_instance}/')
        df_full = pd.read_csv(os.path.join(base_dir, f'{system_name_in_csv_file}-{model_instance}-use_mt_structure=0.csv'))
        out_file = f'pretty_scatterplot-{system_name_in_csv_file}-{model_instance}.png'
        title = f'{system_name_in_csv_file}\n{model_instance}'
        prediction_column = 'substitution_matrix_score'
        ylabel = r'Substitution Matrix Score'
        color = 'tab:red'
    elif model == 'blosum62':
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
        first_wt_mask = np.zeros(df_full.shape[0], dtype=bool)
        first_wt_mask[np.arange(df_full.shape[0])[is_wt_mask][0]] = True
        mask = np.logical_or(np.logical_and(is_not_nan_mask, ~is_wt_mask), first_wt_mask)
    else:
        print('Warning: no "is_wt" column in csv file, assuming there are no wildtype duplicates.', file=sys.stderr)
        mask = is_not_nan_mask

    predictions = df_full[prediction_column][mask]
    targets = df_full[target_column][mask]

    ala_mask = np.logical_or(df_full['scan_type'] == 'alanine', df_full['scan_type'] == 'PIK3CA')
    gly_mask = np.logical_or(df_full['scan_type'] == 'glycine', df_full['scan_type'] == 'PIK3CA')

    pr_all, _ = pearsonr(targets, predictions)
    sr_all, _ = spearmanr(targets, predictions)
    pr_ala, _ = pearsonr(targets[ala_mask], predictions[ala_mask])
    sr_ala, _ = spearmanr(targets[ala_mask], predictions[ala_mask])
    pr_gly, _ = pearsonr(targets[gly_mask], predictions[gly_mask])
    sr_gly, _ = spearmanr(targets[gly_mask], predictions[gly_mask])
    num = targets.shape[0]

    plt.figure(figsize=(3.5, 3.5))

    if args.show_wt_lines in {'both', 'target'}:
        plt.axvline(x=args.wt_value_target, c='k', alpha=1.0, ls='--')

    if args.show_wt_lines in {'both', 'pred'}:
        plt.axhline(y=args.wt_value_pred, c='k', alpha=1.0, ls='--')

    plt.scatter(targets[ala_mask], predictions[ala_mask], c=color, marker='x', alpha=0.5, label='Ala')
    plt.scatter(targets[gly_mask], predictions[gly_mask], c=color, marker='o', alpha=0.5, label='Gly')

    plt.xlabel(args.target_column, fontsize=12)
    plt.ylabel(ylabel, fontsize=12)
    plt.title('Pr. All=%.2f Ala=%.2f Gly=%.2f\nSr. All=%.2f Ala=%.2f Gly=%.2f' % (pr_all, pr_ala, pr_gly, sr_all, sr_ala, sr_gly), fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(base_dir, out_file), dpi=300)
    plt.savefig(os.path.join(base_dir, out_file.replace('.png', '.pdf')), dpi=300)
    plt.close()


