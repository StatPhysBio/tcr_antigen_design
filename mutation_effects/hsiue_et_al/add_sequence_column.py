
import pandas as pd

from constants import WT_SEQ

df = pd.read_csv('hsiue_et_al_H2_sat_mut.csv')

df['sequence'] = [WT_SEQ[:int(mut[1:-1])-1] + mut[-1] + WT_SEQ[int(mut[1:-1]):] for mut in df['mutant']]

df.to_csv('hsiue_et_al_H2_sat_mut.csv', index=False)
