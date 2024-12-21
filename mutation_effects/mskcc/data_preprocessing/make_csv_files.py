
'''

NOTE: there will be repetitions. I include all wildtype entries (e.g. H1H) so as to make it easy to compute the heatmap.
To instead compute a correlation against all values, all wiltype measurements must be collapsed. I will introduce an extra column "is_wt" to more easily indicate that a row is a wiltype row.

'''

import numpy as np
import pandas as pd
# Keeping pandas from truncating long strings
pd.set_option('display.max_colwidth', None)

import json

def reliability_mask(df, mutants):
    is_reliable = []
    columns = ['0.01', '1.0', '100.0']
    for mutant in mutants:
        sub_df = df.loc[df['peptide'] == mutant]
        assert len(sub_df) == 1

        criterion_1 = np.sum(np.isnan(sub_df[columns].values)) == 0
        criterion_2 = np.nanmax(sub_df[columns].values) > 0.06

        is_reliable.append(criterion_1 and criterion_2)
    
    return np.array(is_reliable)



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
wt_peptides = ['NLVPMVATV', 'NLVPMVATV', 'NLVPMVATV', 'IMDQVPFSV', 'IMDQVPFSV', 'IMDQVPFSV', 'GRLKALCQR']

for wt_pdb, pep_chain, pep_resnum_start, tcr, pep_start, wt_pep in zip(wt_pdbs, pep_chains, pep_resnum_start_list, tcrs, pep_start_list, wt_peptides):

    columns = ['mutant', 'is_wt', '-log10(EC50)', 'wt_pdb', 'mt_pdb', 'mutant_chain']

    df = pd.DataFrame(columns=columns)

    with open(f'{pep_start}_{tcr}_ec50s.json', 'r') as f:
        mutant_to_ec50 = json.load(f)

    mutants = []
    is_wt = []
    ec50s = []
    sequences = []
    for mutant in mutant_to_ec50:
        wt_aa = mutant[0]
        resnum = int(mutant[1:-1])
        mt_aa = mutant[-1]

        mutant_in_pdb = f'{wt_aa}{resnum - 1 + pep_resnum_start}{mt_aa}'

        mutants.append(mutant_in_pdb)
        is_wt.append(wt_aa == mt_aa)
        ec50s.append(mutant_to_ec50[mutant])
        sequences.append(f'{wt_pep[:resnum - 1]}{mt_aa}{wt_pep[resnum:]}')

    df['mutant'] = mutants
    df['is_wt'] = is_wt
    df['-log10(EC50)'] = -np.array(ec50s)

    df['wt_seq'] = np.full(len(mutants), wt_pep)

    df['sequence'] = sequences

    # make the reliability mask!
    df['is_reliable'] = reliability_mask(pd.read_csv(f'{pep_start}_{tcr}_reactivity_curves.tsv', sep='\t'), mutants)


    df['wt_pdb'] = np.full(len(mutants), wt_pdb)
    # df['mt_pdb'] = np.full(len(mutants), np.nan)

    df['mutant_chain'] = np.full(len(mutants), pep_chain)

    df.to_csv(f'mskcc_{tcr}_ec50_sat_mut_af3__abs_limit_6.csv', index=False)





