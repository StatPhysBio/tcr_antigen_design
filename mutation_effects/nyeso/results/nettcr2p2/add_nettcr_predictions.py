
import os
import pandas as pd

df = pd.read_csv(f'../../nyeso_peptide_kd_closest.csv')

df_predictions = pd.read_csv('nettcr_predictions_nyeso.csv')
df_predictions = df_predictions.rename(columns={'peptide': 'sequence', 'prediction': 'nettcr_score'})

# Merge the two dataframes on the 'sequence' column
df = pd.merge(df, df_predictions, on='sequence', how='left')

# remove duplicates based upon the "mutant" column
df = df.drop_duplicates(subset=['mutants'])

df.to_csv('nyeso_peptide_kd_closest-nettcr2p2.csv', index=False)

