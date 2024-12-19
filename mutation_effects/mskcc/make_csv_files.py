
'''

NOTE: there will be repetitions. I include all wildtype entries (e.g. H1H) so as to make it easy to compute the heatmap.
To instead compute a correlation against all values, all wiltype measurements must be collapsed. I will introduce an extra column "is_wt" to more easily indicate that a row is a wiltype row.

'''

import numpy as np
import pandas as pd
# Keeping pandas from truncating long strings
pd.set_option('display.max_colwidth', None)

import json

wt_pdbs = [
    '5d2n-filtered',
    'af3_tcr2',
    'af3_tcr3',
    'af3_tcr4',
    'af3_tcr5',
    'af3_tcr6',
    'af3_tcr7'
]
pep_chains = ['I', 'B', 'B', 'B', 'B', 'B', 'B']
pep_resnum_start_list = [1, 1, 1, 1, 1, 1, 1]
tcrs = ['tcr1', 'tcr2', 'tcr3', 'tcr4', 'tcr5', 'tcr6', 'tcr7']
pep_start_list = ['NLV', 'NLV', 'NLV', 'IMD', 'IMD', 'IMD', 'GRL']

for wt_pdb, pep_chain, pep_resnum_start, tcr, pep_start in zip(wt_pdbs, pep_chains, pep_resnum_start_list, tcrs, pep_start_list):

    columns = ['mutant', 'is_wt', '-log10(EC50)', 'wt_pdb', 'mt_pdb', 'mutant_chain']

    df = pd.DataFrame(columns=columns)

    with open(f'{pep_start}_{tcr}_ec50s.json', 'r') as f:
        mutant_to_ec50 = json.load(f)

    mutants = []
    is_wt = []
    ec50s = []
    for mutant in mutant_to_ec50:
        wt_aa = mutant[0]
        resnum = int(mutant[1:-1])
        mt_aa = mutant[-1]

        mutant_in_pdb = f'{wt_aa}{resnum - 1 + pep_resnum_start}{mt_aa}'

        mutants.append(mutant_in_pdb)
        is_wt.append(wt_aa == mt_aa)
        ec50s.append(mutant_to_ec50[mutant])

    df['mutant'] = mutants
    df['is_wt'] = is_wt
    df['-log10(EC50)'] = -np.array(ec50s)


    df['wt_pdb'] = np.full(len(mutants), wt_pdb)
    # df['mt_pdb'] = np.full(len(mutants), np.nan)

    df['mutant_chain'] = np.full(len(mutants), pep_chain)

    df.to_csv(f'luksza_et_al_{tcr}_ec50_sat_mut.csv', index=False)


# # filtered version, excluding the values that are at the upper bound of EC50 which I think are just wrong? Or they mean no binding so in reality they are meaningful? anyway, just doing it to see how it is, it costs me very little
# rows_to_drop = df[df['-EC50'] == -15.549308757970227].index
# df_filtered = df.drop(rows_to_drop, axis=0)
# df_filtered.to_csv('luksza_et_al_tcr1_ec50_sat_mut_FILTERED_HIGH_EC50s.csv', index=False)


