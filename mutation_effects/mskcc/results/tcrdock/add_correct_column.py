
import os
import pandas as pd

tcrs = range(1, 8)

for tcr in tcrs:

    df_true = pd.read_csv(f'../../mskcc_tcr{tcr}_ec50_sat_mut.csv')[['mutant', 'is_wt', '- delta log_ec50_M']]

    df_pred = pd.read_csv(f'_mskcc_tcr{tcr}_ec50_sat_mut_w_pae.tsv', sep='\t')

    # merge
    df = pd.merge(df_true, df_pred, on='mutant', how='inner')

    df.to_csv(f'mskcc_tcr{tcr}_ec50_sat_mut_w_pae.tsv', sep='\t', index=False)



    

