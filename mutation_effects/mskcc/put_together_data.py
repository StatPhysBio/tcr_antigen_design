
import numpy as np
import pandas as pd


wildtypes = {
    '1': 'NLVPMVATV',
    '2': 'NLVPMVATV',
    '3': 'NLVPMVATV',
    '4': 'IMDQVPFSV',
    '5': 'IMDQVPFSV',
    '6': 'IMDQVPFSV',
    '7': 'GRLKALCQR',
}

offsets = {
    '1': 0,
    '2': 375,
    '3': 375,
    '4': 375,
    '5': 375,
    '6': 375,
    '7': 375,
}

dfs = []

for tcr_num in wildtypes:
    df = pd.read_csv(f'mskcc_tcr{tcr_num}_ec50_sat_mut.csv')
    wt_seq = wildtypes[tcr_num]
    offset = offsets[tcr_num]
    resnums = [int(mut[1:-1])-offset for mut in df['mutant']]
    mt_aas = [mut[-1] for mut in df['mutant']]
    df['peptide'] = [wt_seq[:resnum-1]+mt_aa+wt_seq[resnum:] for resnum, mt_aa in zip(resnums, mt_aas)]
    df['system'] = [f'TCR{tcr_num}' for _ in range(len(df))]
    df['study'] = ['luksza_et_al' for _ in range(len(df))]

    dfs.append(df)

df = pd.concat(dfs)

# only keep: peptide, - delta log_ec50_M, system, study
df = df[['peptide', '- delta log_ec50_M', 'system', 'study']]

df.to_csv('peptide_ec50_luksza_et_al.csv', index=False)
