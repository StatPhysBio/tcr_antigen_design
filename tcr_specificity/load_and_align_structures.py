# PyMOL Script: Load and Align Structures by Chain A
# Usage: Modify the 'structure_paths' list to include the paths to your PDB files.

def load_and_align_structures(structure_paths, chain="A"):
    """
    Load multiple structures and align them based on a specified chain.

    Parameters:
        structure_paths (list): List of file paths to the PDB structures.
        reference_name (str): Name for the first loaded structure, used as reference.
        chain (str): Chain identifier to use for alignment (default: "A").
    """

    import os
    from pymol import cmd

    if not structure_paths:
        print("No structure paths provided. Please update the 'structure_paths' list.")
        return

    # Load the first structure and set it as the reference
    reference_name = os.path.basename(structure_paths[0]).split(".")[0]
    cmd.load(structure_paths[0], reference_name)
    cmd.select("reference_chain", f"{reference_name} and chain {chain}")

    # Loop through the remaining structures
    for i, path in enumerate(structure_paths[1:], start=1):
        obj_name = os.path.basename(path).split(".")[0]
        try:
            cmd.load(path, obj_name)
            cmd.select(f"{obj_name}_chain", f"{obj_name} and chain {chain}")
            cmd.align(f"{obj_name}_chain", "reference_chain")
        except:
            print(f"Failed to load and align {path}")

    # Cleanup temporary selections
    cmd.delete("reference_chain")
    for i in range(1, len(structure_paths)):
        cmd.delete(f"structure_{i}_chain")


allele = 'B0702'
pep_length = 9


# import os
# import numpy as np
# import pandas as pd

# df = pd.read_csv('pmhc_class_1_crystal_structures.csv')
# pdbdir = './pdbs/pmhc'

# df_subset = df.loc[np.logical_and(df['mhc_allele_final'] == allele, df['peptide_length'] == pep_length)]

# structure_paths = [os.path.join(pdbdir, pdb+'.pdb') for pdb in df_subset['pdbid']]

# with open(f'structure_paths_{allele}_{pep_length}.txt', 'w+') as f:
#     for path in structure_paths:
#         f.write(path+'\n')


structure_paths = []
with open(f'structure_paths_{allele}_{pep_length}.txt', 'r') as f:
    for line in f:
        structure_paths.append(line.strip())

load_and_align_structures(structure_paths)

