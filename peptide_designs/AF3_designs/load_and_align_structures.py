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


import os

## We assume that chain A is always the MHC chain!

# start with the wildtypes and align everything to them

# ## nyeso
# structure_paths = ['/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/nyeso/pdbs/2bnq.pdb',
#                    '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/nyeso/pdbs/2bnr.pdb']
# for zipfile in os.listdir('nyeso'):
#     if zipfile.startswith('fold') and zipfile.endswith('.zip'):
#         # unzip
#         folder = zipfile[:-4]
#         os.system(f'unzip nyeso/{zipfile} -d nyeso/{folder}')
#         # load structure with best pTM score
#         file = os.path.abspath(f'nyeso/{folder}/{folder}_model_0.cif')
#         structure_paths.append(file)

# ## ebv
# structure_paths = ['/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/hpvg/pdbs/3mv7.pdb.human.MH1.B-35.A.C.DE.pdb',
#                    '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/hpvg_q5/pdbs/4prp.pdb.human.MH1.B-35.A.C.DE.pdb']
# for zipfile in os.listdir('ebv'):
#     if zipfile.startswith('fold') and zipfile.endswith('.zip'):
#         # unzip
#         folder = zipfile[:-4]
#         os.system(f'unzip ebv/{zipfile} -d ebv/{folder}')
#         # load structure with best pTM score
#         file = os.path.abspath(f'ebv/{folder}/{folder}_model_0.cif')
#         structure_paths.append(file)

## magea3
structure_paths = ['/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/magea3/pdbs/5brz.pdb.human.MH1.A-01.A.C.DE.pdb',
                   '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/titin/pdbs/5bs0.pdb.human.MH1.A-01.A.C.DE.pdb']
for zipfile in os.listdir('magea3'):
    if zipfile.startswith('fold') and zipfile.endswith('.zip'):
        # unzip
        folder = zipfile[:-4]
        os.system(f'unzip magea3/{zipfile} -d magea3/{folder}')
        # load structure with best pTM score
        file = os.path.abspath(f'magea3/{folder}/{folder}_model_0.cif')
        structure_paths.append(file)


# Call the function
load_and_align_structures(structure_paths)