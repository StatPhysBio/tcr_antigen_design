
import os
from glob import glob

input_pdbdir = '/gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/all_structures/pdb/ternary'
output_pdbdir = './pdbs/tcr_pmhc'
os.makedirs(output_pdbdir, exist_ok=True)

filenames = os.listdir('/gscratch/spe/gvisan01/tcr_pmhc/tcr_specificity/pwm_csv_files/mhc_motifs/')

pdbids = [filename[:4] for filename in filenames]

for pdbid in pdbids:

    candidate_pdbfiles = glob(os.path.join(input_pdbdir, f'{pdbid}*.pdb'))
    if len(candidate_pdbfiles) == 0:
        print(f'No PDB file found for {pdbid}')
        continue
    elif len(candidate_pdbfiles) > 1:
        print(f'Multiple PDB files found for {pdbid}: {candidate_pdbfiles}')
        continue
    else:
        input_pdbfile = candidate_pdbfiles[0]
    
    output_pdbfile = os.path.join(output_pdbdir, f'{pdbid}.pdb')
    os.system(f'cp {input_pdbfile} {output_pdbfile}')


