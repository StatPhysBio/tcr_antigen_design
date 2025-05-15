
import os
import numpy as np
import pandas as pd


df = pd.read_csv('pmhc_class_1_crystal_structures.csv')


def make_pretty_allele(allele):
    if '-' in allele: # mouse
        pretty_allele = allele
    else:
        pretty_allele = allele[0] + '*' + allele[1:3] + ':' + allele[3:]
    return pretty_allele


allele_and_length_list = df[['mhc_allele_final', 'peptide_length']].values.tolist()
allele_and_length_list = [(x[0], x[1]) for x in allele_and_length_list]
allele_and_length_list = list(set(allele_and_length_list))
allele_and_length_list = sorted(allele_and_length_list, key=lambda x: (x[0], x[1]))

df_out = pd.DataFrame(columns=['MHC allele', 'Peptide Length', 'PDBid'])
for allele_and_length in allele_and_length_list:
    allele, length = allele_and_length
    pretty_allele_name = make_pretty_allele(allele)

    df_curr = df[(df['mhc_allele_final'] == allele) & (df['peptide_length'] == length)]
    pdbids = [x.upper() for x in df_curr['pdbid'].values.tolist()]

    if len(pdbids) >= 2:
        df_out = df_out.append({'MHC allele': pretty_allele_name,
                                'Peptide Length': length,
                                'PDBid': '; '.join(pdbids)}, ignore_index=True)

df_out = df_out.sort_values(by=['MHC allele', 'Peptide Length'])
df_out = df_out.reset_index(drop=True)

# Save the dataframe to a CSV file
df_out.to_csv('pmhc_class_1_crystal_structures_for_paper.csv', index=False)

