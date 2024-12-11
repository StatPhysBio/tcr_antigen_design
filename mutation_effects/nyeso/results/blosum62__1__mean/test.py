
import pandas as pd

from scipy.stats import pearsonr, spearmanr

df = pd.read_csv('nyeso_peptide_kd_2bnq-blosum62__1__mean-use_mt_structure=0__manual.csv')

print('Pearson correlation:', pearsonr(df['-log10(Kd)'], df['substitution_matrix_score']))
print('Spearman correlation:', spearmanr(df['-log10(Kd)'], df['substitution_matrix_score']))
