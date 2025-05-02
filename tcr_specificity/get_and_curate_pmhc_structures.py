
import os, sys
import gzip
import numpy as np
import pandas as pd
from tqdm import tqdm

import argparse

def parse_pdb_chain_lengths_and_resnums(pdb_file_path):
    """
    Parses a PDB file and returns:
    - a dictionary mapping chain IDs to the number of residues in each chain
    - a dictionary mapping chain IDs to a list of residue sequence numbers (resnums), excluding water (HOH)

    Args:
        pdb_file_path (str): Path to the PDB file.

    Returns:
        tuple: 
            - chain_lengths (dict): keys are chain IDs, values are the number of residues
            - chain_resnums (dict): keys are chain IDs, values are lists of residue sequence numbers (as integers)
    """
    chain_lengths = {}
    chain_resnums = {}
    seen_residues = set()

    with open(pdb_file_path, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_name = line[17:20].strip()
                if res_name == "HOH":
                    continue  # Skip water

                chain_id = line[21]
                res_seq = line[22:26].strip()
                try:
                    res_seq_int = int(res_seq)
                except ValueError:
                    continue  # Skip invalid residue numbers

                res_uid = (chain_id, res_seq_int)

                if res_uid not in seen_residues:
                    seen_residues.add(res_uid)
                    chain_lengths[chain_id] = chain_lengths.get(chain_id, 0) + 1
                    chain_resnums.setdefault(chain_id, []).append(res_seq_int)

    return chain_lengths, chain_resnums



if __name__ == '__main__':


    pdbdir = '/gscratch/spe/gvisan01/tcr_pmhc/tcr_specificity/pdbs/pmhc'
    os.makedirs(pdbdir, exist_ok=True)

    df = pd.read_csv('TCR3d_data.csv')
    df['PDB ID'] = df['PDB ID'].apply(lambda x: x.lower())
    df['Species'] = df['Species'].apply(lambda x: x.lower() if isinstance(x, str) else x)
    df = df.rename(columns={'PDB ID': 'pdbid',
                            'Species': 'organism',
                            'MHC allele': 'mhc_allele',
                            'Peptide*': 'wt_peptide',
                            'Bound to TCR': 'bound_to_tcr',
                            'Release date': 'relese_date',
                            'Pubmed': 'pubmed',
                            'Resolution': 'resolution'})
    
    # remove rows where mhc_allele is NaN
    df = df.dropna(subset=['mhc_allele'])

    # remove rows corresponding to structures that originally had the TCR in them
    df = df[df['bound_to_tcr'] == 'unbound']


    # download the first biological assembly of each pdb
    # first, make the csvfile, then run the provided bash script
    # temp_filename = 'pdbs_to_download.txt'
    # with open(temp_filename, 'w+') as f:
    #     f.write(','.join(df['pdbid'].values.tolist()))
    # os.system(f'bash batch_download_pdbs.sh -f {temp_filename} -o {pdbdir} -a')

    # # unzip the pdb files
    # for gz_file in os.listdir(pdbdir):
    #     if gz_file.endswith('.gz'):
    #         with gzip.open(os.path.join(pdbdir, gz_file), 'rt') as f_in:
    #             with open(os.path.join(pdbdir, gz_file.strip('.gz')), 'w+') as f_out:
    #                 f_out.write(f_in.read())
        
    #         # rename from .pdb1 to .pdb
    #         orig_filepath = os.path.join(pdbdir, gz_file.strip('.gz'))
    #         new_filepath = orig_filepath.replace('.pdb1', '.pdb')
    #         os.system(f'mv {orig_filepath} {new_filepath}')

    #         # delete .gz file
    #         os.system(f'rm {os.path.join(pdbdir, gz_file)}')


    # figure out which chain is the peptide chain - same length as the wt peptide in the csv file
    # I should also match the amino-acid sequence to be very extra sure, but I really doubt there are
    # other chains of the same exact length as the peptide that could be confused for the peptide
    print('Identifying peptide chain...', flush=True)
    peptide_chain = []
    peptide_resnums = []
    for i, row in tqdm(df.iterrows(), total=len(df)):
        pdb = row['pdbid']
        expected_peptide_length = len(row['wt_peptide'])

        chain_to_length, chain_to_resnums = parse_pdb_chain_lengths_and_resnums(os.path.join(pdbdir, pdb+'.pdb'))

        for chain in chain_to_length:
            assert chain_to_length[chain] == len(chain_to_resnums[chain])
            if chain_to_length[chain] == expected_peptide_length:
                peptide_chain.append(chain)
                peptide_resnums.append('|'.join(map(str, chain_to_resnums[chain])))
                break
        else:
            print(f'Warning: no chain found that matches the expected length of the peptide for pdb {pdb}.')
            # these pdbs do have the peptide, it's just that it's listed as being on the same chain as the MHC
            # I don't know of any way to do this other than manually on PyMol, which is not worth it for the purposes
            # of this experiment
            peptide_chain.append(np.nan)
            peptide_resnums.append('')

    df['peptide_chain'] = peptide_chain
    df['peptide_resnums'] = peptide_resnums

    # drop rows that are nan in peptide chain
    df = df.dropna(subset=['peptide_chain'])

    df.to_csv('pmhc_class_1_crystal_structures.csv', index=None)


