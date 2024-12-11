
import pandas as pd

from scipy.stats import pearsonr, spearmanr

df = pd.read_csv('mart_peptide_kd_resp_struc_with_relaxation-hermes_py_000-use_mt_structure=0__manual.csv')

print('Pearson correlation:', pearsonr(df['-log10(Kd)'], df['pnE']))
print('Spearman correlation:', spearmanr(df['-log10(Kd)'], df['pnE']))


