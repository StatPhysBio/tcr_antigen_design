
# just for sanity check purposes, which passed

import os
import numpy as np
import pandas as pd

from constants import PDB_TO_WT_SEQ, PDB_TO_PEP_CHAIN, PDBS

df_pred_full = pd.read_csv('hermes_py_000_output.csv')

for pdb in PDBS:

    df_pred_full_pdb = df_pred_full[df_pred_full['pdb'] == pdb]

    df_script = pd.read_csv(f'./results/hermes_py_000/tax_peptide_kd_{pdb}-hermes_py_000-use_mt_structure=0.csv')

    full_preds = []
    print()
    for i, row in df_script.iterrows():
        mutants = row['mutants'].split('|')
        preds = []
        for mutant in mutants:
            wt_aa = mutant[0]
            mt_aa = mutant[-1]
            pos = int(mutant[1:-1])

            wt_logproba = df_pred_full_pdb[df_pred_full_pdb['resnum'] == pos][f'logproba_{wt_aa}'].values[0]
            mt_logproba = df_pred_full_pdb[df_pred_full_pdb['resnum'] == pos][f'logproba_{mt_aa}'].values[0]

            # print(mt_logproba, wt_logproba)

            preds.append(mt_logproba - wt_logproba)
        
        mean_pred = np.mean(preds)
        full_preds.append(mean_pred)

        script_pred = row['log_proba_mt__minus__log_proba_wt']

        print('%s\t%.4f\t%.4f\t%s' % (pdb, script_pred, mean_pred, mutants))
    




