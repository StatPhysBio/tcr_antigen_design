
import os
import gzip, pickle
import numpy as np
import pandas as pd


with gzip.open('/gscratch/spe/gvisan01/peptide_mhc/mhc_motif_atlas/mhc_motif_pwms.pkl.gz', 'rb') as f:
    mhc_motif_pwms = pickle.load(f)

print(mhc_motif_pwms['class_I'].keys())

df = pd.read_csv('pmhc_class_1_crystal_structures.csv')

def fix_allele(allele):
    
    if isinstance(allele, str):
        if '*' in allele:
            return allele[:1] + ''.join(allele[2:].split(':')[:2])
        else:
            return allele[:2] + '-' + allele[2:]
    else:
        return allele

df['mhc_allele_final'] = df['mhc_allele_longer'].apply(fix_allele)



