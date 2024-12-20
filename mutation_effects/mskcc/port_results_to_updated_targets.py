
import os
import numpy as np
import pandas as pd

identifier_template = 'mskcc_tcr{tcr}_ec50_sat_mut_af3'

tcrs = [1, 2, 3, 4, 5, 6, 7]


## 1) Hermes models

columns_to_drop = ['Unnamed: 0', '- delta log_ec50_M']
columns_to_merge_over = ['mutant', 'sequence', 'wt_seq', 'wt_pdb', 'is_wt', 'mutant_chain']
columns_to_merge_over_backup = ['mutant', 'sequence', 'wt_pdb', 'is_wt', 'mutant_chain']

models = ['hermes_py_000', 'hermes_py_050']

for model in models:

    for tcr in tcrs:

        identifier = identifier_template.format(tcr=tcr)

        for relaxation_str in ['', '_with_relaxation']:

            df_new_targets = pd.read_csv(f'{identifier}.csv')

            try:
                df_old_pred = pd.read_csv(f'results_old/{model}/{identifier}{relaxation_str}-{model}-use_mt_structure=0.csv')
            except FileNotFoundError:
                continue

            df_old_pred = df_old_pred.drop(columns=columns_to_drop)

            try:
                df_merged = pd.merge(df_new_targets, df_old_pred, on=columns_to_merge_over, how='inner')
            except:
                df_merged = pd.merge(df_new_targets, df_old_pred, on=columns_to_merge_over_backup, how='inner')
            
            os.makedirs(f'results/{model}', exist_ok=True)
            df_merged.to_csv(f'results/{model}/{identifier}{relaxation_str}-{model}-use_mt_structure=0.csv', index=False)


## 2) proteinmpnn models

columns_to_drop = ['Unnamed: 0', '- delta log_ec50_M']
columns_to_merge_over = ['mutant', 'sequence', 'wt_seq', 'wt_pdb', 'is_wt', 'mutant_chain']
columns_to_merge_over_backup = ['mutant', 'sequence', 'wt_pdb', 'is_wt', 'mutant_chain']

models = ['proteinmpnn_v_48_002', 'proteinmpnn_v_48_030']

for model in models:

    for tcr in tcrs:

        identifier = identifier_template.format(tcr=tcr)

        df_new_targets = pd.read_csv(f'{identifier}.csv')

        try:
            df_old_pred = pd.read_csv(f'results_old/{model}/zero_shot_predictions/{identifier}-num_seq_per_target=10-use_mt_structure=0.csv')
        except FileNotFoundError:
            continue

        df_old_pred = df_old_pred.drop(columns=columns_to_drop)

        try:
            df_merged = pd.merge(df_new_targets, df_old_pred, on=columns_to_merge_over, how='inner')
        except:
            df_merged = pd.merge(df_new_targets, df_old_pred, on=columns_to_merge_over_backup, how='inner')
        
        os.makedirs(f'results/{model}/zero_shot_predictions', exist_ok=True)
        df_merged.to_csv(f'results/{model}/zero_shot_predictions/{identifier}-num_seq_per_target=10-use_mt_structure=0.csv', index=False)


## 3) blosum62

columns_to_drop = ['Unnamed: 0', '- delta log_ec50_M']
columns_to_merge_over = ['mutant', 'sequence', 'wt_seq', 'wt_pdb', 'is_wt', 'mutant_chain']
columns_to_merge_over_backup = ['mutant', 'sequence', 'wt_pdb', 'is_wt', 'mutant_chain']

models = ['blosum62']

for model in models:

    for tcr in tcrs:

        identifier = identifier_template.format(tcr=tcr)

        df_new_targets = pd.read_csv(f'{identifier}.csv')

        try:
            df_old_pred = pd.read_csv(f'results_old/{model}/{identifier}-{model}-use_mt_structure=0.csv')
        except FileNotFoundError:
            continue

        df_old_pred = df_old_pred.drop(columns=columns_to_drop)

        try:
            df_merged = pd.merge(df_new_targets, df_old_pred, on=columns_to_merge_over, how='inner')
        except:
            df_merged = pd.merge(df_new_targets, df_old_pred, on=columns_to_merge_over_backup, how='inner')
        
        os.makedirs(f'results/{model}', exist_ok=True)
        df_merged.to_csv(f'results/{model}/{identifier}-{model}-use_mt_structure=0.csv', index=False)


## 4) TCRdock

columns_to_drop = ['- delta log_ec50_M', '-EC50', 'mutant', 'mutant_chain'] # mutant and mutant_chain were with tcrdock numbering so wrong
columns_to_merge_over = ['peptide', 'is_wt']

models = ['tcrdock', 'tcrdock_no_nearby_templates']

for model in models:

    for tcr in tcrs:

        identifier = identifier_template.format(tcr=tcr)

        df_new_targets = pd.read_csv(f'{identifier}.csv')
        df_new_targets = df_new_targets.rename(columns={'sequence': 'peptide'})

        try:
            if model == 'tcrdock':
                df_old_pred = pd.read_csv(f'results_old/{model}/{identifier}_w_pae.tsv', sep='\t')
            else:
                df_old_pred = pd.read_csv(f'results_old/{model}/{identifier}_no_nearby_templates_w_pae.tsv', sep='\t')
        except FileNotFoundError:
            continue
        
        df_old_pred = df_old_pred.drop(columns=columns_to_drop)

        df_merged = pd.merge(df_new_targets, df_old_pred, on=columns_to_merge_over, how='inner')

        # only keep one sintance of 'is_wt'
        df_merged = df_merged.drop_duplicates(subset=['peptide', 'is_wt'])

        os.makedirs(f'results/{model}', exist_ok=True)
        df_merged.to_csv(f'results/{model}/{identifier}_no_nearby_templates_w_pae.tsv', sep='\t', index=False)
