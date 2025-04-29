
import os
from glob import glob
import numpy as np
import pandas as pd

for tcr in range(1, 7+1):

    files = glob(f'/Users/gianmarcovisani/Downloads/output_mskcc_tcr{tcr}*.tsv')

    df_list = []

    for file in files:

        df_ = pd.read_csv(file, sep='\t')

        df_list.append(df_)

    df_results = pd.concat(df_list, ignore_index=True)
    df_results.rename(columns={'tapir_antigen': 'sequence',
                            'tapir_scores': 'tapir_score'}, inplace=True)

    df_data = pd.read_csv(f'../../mskcc_tcr{tcr}_ec50_sat_mut_af3.csv')

    df_results = pd.merge(df_data, df_results, on=['sequence'], how='left')

    # remove duplicate rows
    df_results = df_results.loc[~df_results.duplicated(subset=['mutant'], keep='first')]

    df_results.to_csv(f'mskcc_tcr{tcr}_ec50_sat_mut_af3-tapir.csv')

