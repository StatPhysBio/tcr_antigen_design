
import os
import pandas as pd


for tcr in range(1, 8):

    df = pd.read_csv(f'../../mskcc_tcr{tcr}_ec50_sat_mut_af3.csv')

    df_predictions = pd.read_csv(f'nettcr_predictions_tcr{tcr}.csv')
    df_predictions = df_predictions.rename(columns={'peptide': 'sequence', 'prediction': 'nettcr_score'})

    # Merge the two dataframes on the 'sequence' column
    df = pd.merge(df, df_predictions, on='sequence', how='left')

    # remove duplicates based upon the "mutant" column
    df = df.drop_duplicates(subset=['mutant'])

    df.to_csv(f'mskcc_tcr{tcr}_ec50_sat_mut_af3-nettcr2p2.csv', index=False)

