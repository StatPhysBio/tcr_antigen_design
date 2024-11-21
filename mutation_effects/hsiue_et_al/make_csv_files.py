
'''

NOTE: there will be repetitions. I include all wildtype entries (e.g. H1H) so as to make it easy to compute the heatmap.
To instead compute a correlation against all values, all wiltype measurements must be collapsed. I will introduce an extra column "is_wt" to more easily indicate that a row is a wiltype row.

'''

import numpy as np
import pandas as pd
# Keeping pandas from truncating long strings
pd.set_option('display.max_colwidth', None)

ifn_gamma_table = pd.read_csv('hsiue_et_al_H2_figure_4D.csv', index_col=0) # some measurement of TCR-peptide avidity (aka activity)
mfi_table = pd.read_csv('hsiue_et_al_H2_figure_S11A.csv', index_col=0) # some measurement of peptide loading efficiency

wt_pep = 'HMTEVVRHC'
pep_len = len(wt_pep)
pep_chain = 'F'
aa_order = 'AILMFVPGWRHKDENQSTYC'
wt_pdb = '6w51-filtered'


columns = ['mutant', 'is_wt', 'IFN_gamma (pg/ml)', 'wt_pdb', 'mt_pdb', 'mutant_chain']

df = pd.DataFrame(columns=columns)

mutants = []
is_wt = []
ifn_gamma = []
mfi = []
for i, wt_aa in enumerate(wt_pep):
    for mt_aa in aa_order:
        mutants.append(wt_aa + str(i+1) + mt_aa)
        if wt_aa == mt_aa:
            is_wt.append(True)
        else:
            is_wt.append(False)
        ifn_gamma.append(float(ifn_gamma_table[mt_aa][i+1]))
        mfi.append(float(mfi_table[mt_aa][i+1]))
        

df['mutant'] = mutants
df['is_wt'] = is_wt

df['IFN_gamma (pg/ml)'] = ifn_gamma
df['MFI'] = mfi

df['wt_pdb'] = np.full(len(mutants), wt_pdb)
# df['mt_pdb'] = np.full(len(mutants), np.nan)

df['mutant_chain'] = np.full(len(mutants), pep_chain)

df.to_csv('hsiue_et_al_H2_sat_mut.csv', index=False)










