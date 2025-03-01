
import numpy as np
import pandas as pd

df = pd.read_csv('nyeso_peptide_kd_closest.csv', index_col=False)

print(-np.log10(df['Kd']))
print(df['-log10(Kd)'])
