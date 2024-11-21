
import numpy as np
import pandas as pd

from constants import PDB_TO_WT_SEQ, PDB_TO_PEP_CHAIN

def get_mutants(wt_seq, mt_seq):
    assert len(wt_seq) == len(mt_seq)
    return [f'{wt}{i+1}{mt}' for i, (wt, mt) in enumerate(zip(wt_seq, mt_seq)) if wt != mt]


df = pd.read_csv('nyeso_peptide_kd.csv')

columns_to_keep = 'Kd,system,name,sequence,-log10(Kd),study'.split(',')

for pdb, wt_seq in PDB_TO_WT_SEQ.items():
    
    curr_df = df[columns_to_keep].copy()

    mutants = ['|'.join(get_mutants(wt_seq, x)) for x in curr_df['sequence'].values]
    new_mutants = []
    for mutant in mutants:
        if mutant:
            new_mutants.append(mutant)
        else:
            new_mutants.append(f'{wt_seq[0]}1{wt_seq[0]}')
    curr_df['mutants'] = new_mutants

    curr_df['wt_seq'] = wt_seq
    curr_df['pdb'] = pdb
    curr_df['chain'] = PDB_TO_PEP_CHAIN[pdb]

    curr_df.to_csv(f'nyeso_peptide_kd_{pdb}.csv', index=False)


