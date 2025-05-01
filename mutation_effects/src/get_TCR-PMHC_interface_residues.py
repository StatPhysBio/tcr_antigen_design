
import os
import numpy as np
import pandas as pd
import json

from Bio.PDB import PDBParser, is_aa
from Bio.PDB.NeighborSearch import NeighborSearch

import argparse



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--pdbfile', type=str, required=True)
    parser.add_argument('-p', '--peptide_chain', type=str, required=True)
    parser.add_argument('-a', '--tcrA_chain', type=str, required=True)
    parser.add_argument('-b', '--tcrB_chain', type=str, required=True)
    args = parser.parse_args()

    # Load the structure
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('protein', args.pdbfile)

    model = structure[0]

    def get_valid_ca_atoms(chain):
        """Return list of (residue, CA atom) for standard amino acids (excluding water)."""
        return [
            (res, res['CA']) for res in chain
            if is_aa(res, standard=True) and res.get_resname() != 'HOH' and 'CA' in res
        ]

    # Get alpha-carbons from chain P
    ca_atoms_P = get_valid_ca_atoms(model[args.peptide_chain])
    ca_atom_objs_P = [atom for (_, atom) in ca_atoms_P]
    ns = NeighborSearch(ca_atom_objs_P)

    selected_residues_A = set()
    selected_residues_B = set()
    selected_residues_P = set()
    all_residues_P = set(res.get_id()[1] for (res, _) in ca_atoms_P)

    # Check chains A and B
    for chain_id in [args.tcrA_chain, args.tcrB_chain]:
        ca_atoms_chain = get_valid_ca_atoms(model[chain_id])
        for res, ca in ca_atoms_chain:
            nearby = ns.search(ca.coord, 12.0, level='A')
            if nearby:
                if chain_id == args.tcrA_chain:
                    selected_residues_A.add(res.get_id()[1])
                else:
                    selected_residues_B.add(res.get_id()[1])
                # Add the involved P residues (whose CA was found)
                for atom in nearby:
                    resP = atom.get_parent()
                    selected_residues_P.add(resP.get_id()[1])

    # save output as a JSON
    output_dict = {
        'all_peptide_residues': {args.peptide_chain: list(all_residues_P)},
        'interacting_peptide_residues': {args.peptide_chain: list(selected_residues_P)},
        'interacting_TCR_residues': {args.tcrA_chain: list(selected_residues_A),
                                     args.tcrB_chain: list(selected_residues_B)},
    }

    with open(args.pdbfile.replace('.pdb', '__tcr_pep_interactions.json'), 'w') as f:
        json.dump(output_dict, f, indent=4)


