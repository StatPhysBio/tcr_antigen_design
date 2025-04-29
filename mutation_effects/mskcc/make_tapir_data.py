
import os
import numpy as np
import pandas as pd

for tcr in range(1, 8):

    df_in = pd.read_csv(f'results/tcrdock/mskcc_tcr{tcr}_ec50_sat_mut_af3_w_pae.tsv', sep='\t')

    df_out = pd.DataFrame({'alpha_v': [df_in['va'].values[0][:-3]],  # removing the *01
            'alpha_j': [df_in['ja'].values[0][:-3]],  # removing the *01
            'alpha_cdr3': [df_in['cdr3a'].values[0]], 
            'beta_v': [df_in['vb'].values[0][:-3]],  # removing the *01
            'beta_j': [df_in['jb'].values[0][:-3]],  # removing the *01
            'beta_cdr3': [df_in['cdr3b'].values[0]]})

    df_out.to_csv(f'mskcc_tcr{tcr}_tapir.tsv', sep='\t', index=False)




