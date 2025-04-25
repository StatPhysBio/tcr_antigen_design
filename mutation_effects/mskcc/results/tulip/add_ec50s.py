
import os
import sys
import numpy as np
import pandas as pd


for tcr in ['1', '2', '3', '4', '5', '6']:
    # just merge!

    df = pd.read_csv(f'tcr{tcr}_score.csv')
    df.rename(columns={'peptide': 'sequence'}, inplace=True)

    df_with_ec_50 = pd.read_csv(f'../../mskcc_tcr{tcr}_ec50_sat_mut_af3.csv')

    df = df.merge(df_with_ec_50, on='sequence', how='left')

    df.to_csv(f'mskcc_tcr{tcr}_ec50_sat_mut_af3-tulip.csv', index=False)

