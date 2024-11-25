
import numpy as np
import pandas as pd

from constants import PDB_TO_WT_SEQ, PDB_TO_PEP_CHAIN

def get_mutants(wt_seq, mt_seq):
    assert len(wt_seq) == len(mt_seq)
    return [f'{wt}{i+1}{mt}' for i, (wt, mt) in enumerate(zip(wt_seq, mt_seq)) if wt != mt]


df = pd.read_csv('mart_peptide_kd.csv')

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
    curr_df['chain'] = ['|'.join([PDB_TO_PEP_CHAIN[pdb] for _ in range(len(mutant.split('|')))]) for mutant in curr_df['mutants'].values]

    # average the -log10(Kd) values for the same mutant
    curr_df['-log10(Kd)'] = curr_df.groupby('mutants')['-log10(Kd)'].transform('mean')

    # for those rows, make study the concatenation of the respective studies that have been averaged over
    curr_df['study'] = curr_df.groupby('mutants')['study'].transform(lambda x: '|'.join(x))

    # remove duplicate rows (according to mutants)
    curr_df = curr_df.drop_duplicates(subset='mutants')

    # make a new column called "neg_delta_log10_Kd" which is the -log10(Kd) value minus the -log10(Kd) value of the wildtype
    curr_df['neg_delta_log10_Kd'] = curr_df['-log10(Kd)'] - curr_df['-log10(Kd)'][curr_df['sequence'] == wt_seq].values

    curr_df.to_csv(f'mart_peptide_kd_{pdb}.csv', index=False)



df = pd.read_csv('results/tcrdock/mart_peptide_kd_w_pae.tsv', sep='\t')

# average the -log10(Kd) values for the same mutant
df['-log10(Kd)'] = df.groupby('peptide')['-log10(Kd)'].transform('mean')

# remove duplicate rows (according to peptide)
df = df.drop_duplicates(subset='peptide')

df.to_csv('results/tcrdock/mart_peptide_kd_w_pae_filtered.tsv', sep='\t', index=False)




df = pd.read_csv('results/tcrdock_no_nearby_templates/mart_peptide_kd_no_nearby_templates_w_pae.tsv', sep='\t')

# average the -log10(Kd) values for the same mutant
df['-log10(Kd)'] = df.groupby('peptide')['-log10(Kd)'].transform('mean')

# remove duplicate rows (according to peptide)
df = df.drop_duplicates(subset='peptide')

df.to_csv('results/tcrdock_no_nearby_templates/mart_peptide_kd_no_nearby_templates_w_pae_filtered.tsv', sep='\t', index=False)

