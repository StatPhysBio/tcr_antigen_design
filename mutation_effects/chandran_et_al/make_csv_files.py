

import numpy as np
import pandas as pd
# Keeping pandas from truncating long strings
pd.set_option('display.max_colwidth', None)

# # The values below were eye-balled from the figures in the paper. they should be roughly correct
# PRELIM_SCANS_TABLE = {
#     'TCR4': {
#         'alanine': [np.nan, 0.15, 0.20, 0.60, 0, 0, 0, 0, 0],
#         'glycine': [1.05, 0.15, 0, np.nan, np.nan, 0, 0, 0, 0]
#     },
#     'TCR3': {
#         'alanine': [np.nan, 0.80, 0.75, 0.75, 0, 0, 0, 0.25, 0],
#         'glycine': [0.90, 0.35, 0, np.nan, np.nan, 0, 0, 1.00, 0]
#     }
# }

WT_PDB_TABLE = {
    'TCR4': '7l1d-fixed',
    'TCR3': '7rrg-fixed'
}

from pzfx_parser import read_pzfx

raw_data = read_pzfx('Ala and Gly Scans for PIK3CA TCRs.pzfx')

data_keys = {
    'alanine': 'Alanine Scan',
    'glycine': 'Glycine Scan'
}

scan_keys = {
    'alanine': ['P1-x', 'P2A', 'P3A', 'P4A', 'P5A', 'P6A', 'P7A', 'P8A', 'P9A'],
    'glycine': ['P1G', 'P2G', 'P3G', 'P4-x', 'P5-x', 'P6G', 'P7G', 'P8G', 'P9G']
}

tcr_columns = {
    'TCR3': ['TCR3_3', 'TCR3_4', 'TCR3_5'],
    'TCR4': ['TCR4_6', 'TCR4_7', 'TCR4_8']
}

wt_pep = 'ALHGGWTTK'
pep_len = len(wt_pep)
pep_chain = 'C'
orig_wt_mutant = 'L2H'
orig_wt_row = 'WT'


columns = ['mutant', 'TFN_alpha_mt/TFN_alpha_wt', 'wt_pdb', 'mt_pdb', 'mutant_chain', 'TCR', 'scan_type']

df = pd.DataFrame(columns=columns)

all_mutants = []
all_tfn_mt_over_tfn_wt = []
all_wt_pdbs = []
# all_mt_pdbs = []
all_mutant_chains = []
all_tcrs = []
all_scan_types = []


for tcr in ['TCR4', 'TCR3']:

    wt_pdb = WT_PDB_TABLE[tcr]

    for scan_type, mt_aa in zip(['alanine', 'glycine'], ['A', 'G']):

        data_df = raw_data[data_keys[scan_type]]

        mutants = []
        tfn_mt_over_tfn_wt = []
        for i in range(pep_len):
            wt_aa = wt_pep[i]
            resnum = i + 1

            row = data_df.loc[scan_keys[scan_type][i]]

            tfn_ratio = np.mean([np.float32(row[tcr_col]) for tcr_col in tcr_columns[tcr]])

            if np.isnan(tfn_ratio):
                continue
            else:
                mutants.append(f'{wt_aa}{resnum}{mt_aa}')
                tfn_mt_over_tfn_wt.append(tfn_ratio)
        
        if scan_type == 'alanine':
            # add the value for the original PIK3CA (un-mutated), but just do it once
            mutants.append(orig_wt_mutant)
            tfn_mt_over_tfn_wt.append(np.mean([data_df.loc[orig_wt_row][tcr_col] for tcr_col in tcr_columns[tcr]]))
            all_scan_types.extend(list(np.full(len(mutants)-1, scan_type)) + ['PIK3CA'])
        
        all_mutants.extend(mutants)
        all_tfn_mt_over_tfn_wt.extend(tfn_mt_over_tfn_wt)
        all_wt_pdbs.extend(np.full(len(mutants), wt_pdb))
        # all_mt_pdbs.extend(np.full(len(mutants), np.nan))
        all_mutant_chains.extend(np.full(len(mutants), pep_chain))
        all_tcrs.extend(np.full(len(mutants), tcr))
        if scan_type != 'alanine':
            all_scan_types.extend(np.full(len(mutants), scan_type))
        

df['mutant'] = all_mutants
df['TFN_alpha_mt/TFN_alpha_wt'] = all_tfn_mt_over_tfn_wt
df['wt_pdb'] = all_wt_pdbs
# df['mt_pdb'] = all_mt_pdbs
df['mutant_chain'] = all_mutant_chains
df['TCR'] = all_tcrs
df['scan_type'] = all_scan_types
df.to_csv(f'chandran_et_al_peptide_A_and_G_scans.csv', index=False)


