
import os
import pandas as pd

filenames = os.listdir('/gscratch/spe/gvisan01/tcr_pmhc/tcr_specificity/pwm_csv_files/mhc_motifs/')

pdbids = []
wt_peps = []
alleles = []
for filename in filenames:
    filename = filename.split('.')[0]
    pdbid, wt_pep, allele = filename.split('__')

    pdbids.append(pdbid)
    wt_peps.append(wt_pep)
    alleles.append(allele)

df = pd.DataFrame({'pdbid': pdbids, 'wt_peptide': wt_peps, 'mhc_allele': alleles})
df.to_csv('pdbs_allele_df.csv', index=False)
