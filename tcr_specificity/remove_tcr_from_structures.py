
import os
import numpy as np
import pandas as pd

input_pdbdir = './pdbs/tcr_pmhc'
output_pdbdir = './pdbs/pmhc'
os.makedirs(output_pdbdir, exist_ok=True)

# if chain E exists, remove atoms belonging to either chain D or E, as those are always the TCR atoms in these standardized structures,
# if chain E does not exist, remove atoms belonging to chain C or D, as those are always the TCR atoms in these standardized structures

pdbids = []
peptide_chains = []
peptide_resnums = []

for pdbfile in os.listdir(input_pdbdir):
    if pdbfile.endswith('.pdb'):
        input_pdbfile = os.path.join(input_pdbdir, pdbfile)
        output_pdbfile = os.path.join(output_pdbdir, pdbfile)
        
        with open(input_pdbfile, 'r') as f:
            lines = f.readlines()
        
        # check if chain E exists
        chain_e_exists = False
        for line in lines:
            if line.startswith('ATOM') and line[21] == 'E':
                chain_e_exists = True
                break
        
        if chain_e_exists:
            chains_to_remove = ['D', 'E']
            peptide_chain = 'C'
        else:
            chains_to_remove = ['C', 'D']
            peptide_chain = 'B'
        
        curr_resnums = []
        
        
        with open(output_pdbfile, 'w') as f:
            for line in lines:
                if line.startswith('ATOM'):
                    chain_id = line[21]
                    if chain_id not in chains_to_remove:
                        f.write(line)
                else:
                    f.write(line)
                
                # gather the resnums of the peptide chain
                if line.startswith('ATOM') and line[21] == peptide_chain:
                    resnum = line[22:26].strip()
                    if resnum not in curr_resnums:
                        curr_resnums.append(resnum)
                
    
        # remove the TER lines from the end of the file, which are not needed anymore
        with open(output_pdbfile, 'r') as f:
            lines = f.readlines()
        with open(output_pdbfile, 'w') as f:
            done = False
            newlines = []
            for line in lines[::-1]:
                if not line.startswith('TER'):
                    newlines.append(line)
                    done = True
                elif done:
                    newlines.append(line)
            
            for line in newlines[::-1]:
                f.write(line)
        
        pdbids.append(pdbfile[:-4])
        peptide_chains.append(peptide_chain)
        peptide_resnums.append('|'.join(curr_resnums))

curr_df = pd.DataFrame({
    'pdbid': pdbids,
    'peptide_chain': peptide_chains,
    'peptide_resnums': peptide_resnums
})

df_with_other_info = pd.read_csv('pdbs_allele_df.csv')

# merge on the pdbid column
df = pd.merge(curr_df, df_with_other_info, on='pdbid', how='left')

# save the dataframe to a csv file
df.to_csv('pdbs_allele_df_with_peptide_info.csv', index=False)


