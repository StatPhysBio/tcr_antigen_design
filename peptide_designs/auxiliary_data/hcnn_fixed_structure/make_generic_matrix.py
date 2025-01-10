
import os
import numpy as np
import pandas as pd

AMINOACIDS = 'CMFILVWYAGTSNQDEHRKP'


df = pd.read_csv('/gscratch/spe/gvisan01/hermes/experiments/aa_cls_eval/output/hermes_py_000.csv')

matrix_num = np.zeros((20, 20))
matrix_denom = np.zeros((20, 20))

for i, row in df.iterrows():

    aa_wt = row['resname']

    for j, aa_mt in enumerate(AMINOACIDS):
        logit_wt = f'logit_{aa_wt}'
        logit_mt = f'logit_{aa_mt}'
        matrix_num[AMINOACIDS.index(aa_wt), AMINOACIDS.index(aa_mt)] += row[logit_mt] - row[logit_wt]
        matrix_denom[AMINOACIDS.index(aa_wt), AMINOACIDS.index(aa_mt)] += 1

matrix = matrix_num / matrix_denom

matrix_df = pd.DataFrame(matrix, index=list(AMINOACIDS), columns=list(AMINOACIDS))
matrix_df.to_csv(f'hermes_000_on_general_proteins_marginal_matrix.csv')

