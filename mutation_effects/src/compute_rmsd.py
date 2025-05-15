
import os
import numpy as np
from Bio.PDB import PDBParser, is_aa, Superimposer
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.pairwise2 import format_alignment

def align_chains(pdb_file1, chain_id1, resnums1, pdb_file2, chain_id2, resnums2):
    """
    Aligns chain from pdb_file2 to chain from pdb_file1 using CA atoms.

    Returns the superimposer object and aligned structure 2.
    """
    parser = PDBParser(QUIET=True)
    structure1 = parser.get_structure("struct1", pdb_file1)
    structure2 = parser.get_structure("struct2", pdb_file2)

    # Extract chains
    chain1 = structure1[0][chain_id1]
    chain2 = structure2[0][chain_id2]

    # Get CA atoms
    atoms1 = [res['CA'] for res in chain1 if 'CA' in res if res.id[1] in set(resnums1)]
    atoms2 = [res['CA'] for res in chain2 if 'CA' in res if res.id[1] in set(resnums2)]

    # Superimpose
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(structure2.get_atoms())  # Apply transform to structure2

    return super_imposer, structure1, structure2

def align_by_atoms(structure2, atoms1, atoms2):
    """
    Aligns two structures based on given atoms.
    """
    super_imposer = Superimposer()
    super_imposer.set_atoms(atoms1, atoms2)
    super_imposer.apply(structure2.get_atoms())  # Apply transform to structure2
    return super_imposer, structure2


def align_sequences_with_gap_penalty(seq1, seq2, mode='global', match=1, mismatch=-1, gap_open=-10, gap_extend=-0.5):
    """
    Aligns two sequences using Biopython's pairwise2 with gap penalties.

    Parameters:
        seq1, seq2: sequences (strings)
        mode: 'global' or 'local'
        match: score for a match
        mismatch: penalty for a mismatch
        gap_open: penalty for opening a gap
        gap_extend: penalty for extending a gap

    Returns:
        A tuple (aligned_seq1, aligned_seq2, score)
    """
    seq1 = str(seq1)
    seq2 = str(seq2)

    if mode == 'global':
        alignments = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap_open, gap_extend, one_alignment_only=True)
    elif mode == 'local':
        alignments = pairwise2.align.localms(seq1, seq2, match, mismatch, gap_open, gap_extend, one_alignment_only=True)
    else:
        raise ValueError("mode must be 'global' or 'local'")

    aligned_seq1, aligned_seq2, score = alignments[0][:3]
    return aligned_seq1, aligned_seq2, score

def align_sequences_middle_gap_penalty(seq1, seq2, mode='global', match=2, mismatch=-0.1, gap_open=-100, gap_extend=-100):
    """
    Aligns two sequences, penalizing internal gaps more than terminal gaps.

    Parameters:
        seq1, seq2: sequences (strings)
        mode: 'global' or 'local'
        match: match score
        mismatch: mismatch penalty
        gap_open: gap opening penalty
        gap_extend: gap extension penalty

    Returns:
        A tuple (aligned_seq1, aligned_seq2, score)
    """
    seq1 = str(seq1)
    seq2 = str(seq2)

    # Avoid penalizing end gaps
    penalize_ends = (True, True)
    if mode == 'global':
        # Middle gaps are penalized, terminal gaps are not
        alignments = pairwise2.align.globalms(
            seq1, seq2,
            match, mismatch,
            gap_open, gap_extend,
            penalize_end_gaps=(False, False),  # <=== key change
            one_alignment_only=True
        )
    elif mode == 'local':
        alignments = pairwise2.align.localms(
            seq1, seq2,
            match, mismatch,
            gap_open, gap_extend,
            one_alignment_only=True
        )
    else:
        raise ValueError("mode must be 'global' or 'local'")

    aligned_seq1, aligned_seq2, score = alignments[0][:3]
    return aligned_seq1, aligned_seq2, score


def extract_chain_residue_info(pdb_file):
    """
    Extract residue information from a PDB file using Biopython.

    Parameters
    ----------
    pdb_file : str
        Path to the input PDB file.

    Returns
    -------
    dict
        A nested dictionary with the following structure:

        {
            'chain_id': {
                'resnames': [list of residue names],
                'resnums': [list of residue numbers],
                'inscodes': [list of insertion codes]
            },
            ...
        }

        - Residue names are 1-letter codes (e.g., 'A', 'G').
        - Residue numbers are integers.
        - Insertion codes are single characters ('' if absent).

    Notes
    -----
    - Only standard amino acid residues are included.
    - Waters, ligands, and other non-protein residues are excluded.
    - Empty chains (chains with no standard residues) are omitted.
    """

    aa_to_one_letter = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
        'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
        'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
        'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
        'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', pdb_file)

    chain_dict = {}

    for model in structure:
        for chain in model:
            chain_id = chain.id
            resnames = []
            resnums = []
            inscodes = []

            for residue in chain:
                # Skip hetatm and water
                if not is_aa(residue, standard=True):
                    continue

                if residue.get_resname() in aa_to_one_letter:
                    resname = aa_to_one_letter[residue.get_resname()]
                else:
                    resname = 'X'
                resnum = residue.id[1]
                inscode = residue.id[2].strip()  # Insertion code

                resnames.append(resname)
                resnums.append(resnum)
                inscodes.append(inscode)

            # Only store chains with residues
            if resnames:
                chain_dict[chain_id] = {
                    'resnames': resnames,
                    'resnums': resnums,
                    'inscodes': inscodes
                }

    return chain_dict


def build_chain_sequences(chain_dict):
    """
    Build chain sequences with gaps from extracted residue information.

    Parameters
    ----------
    chain_dict : dict
        Output dictionary from `extract_chain_residue_info()`.

    Returns
    -------
    dict
        Dictionary mapping chain IDs to sequences (str),
        where gaps are represented by '-' for missing residues.
    """
    seq_dict = {}

    for chain_id, resdata in chain_dict.items():
        resnames = resdata['resnames']
        resnums = resdata['resnums']
        
        if not resnames:
            seq_dict[chain_id] = ''
            continue
        
        sequence = ''
        # Start from first residue number
        # prev_num = resnums[0] - 1
        
        for resname, resnum in zip(resnames, resnums):
            # gap = resnum - prev_num - 1  # How many residues missing?
            # if gap > 0:
            #     sequence += '-' * gap  # Insert gaps
            sequence += resname
            # prev_num = resnum
        
        seq_dict[chain_id] = sequence
    
    return seq_dict

def get_ca_atoms(structure, chain_id, residue_numbers):
    """
    Given a PDB file, chain ID, and a list of residue numbers,
    return a list of CA atom objects for those residues.

    Parameters:
        structure (PDB.Structure): Biopython Structure object
        chain_id (str): Chain identifier (e.g., 'A')
        residue_numbers (list of int): List of residue numbers to extract

    Returns:
        list of Bio.PDB.Atom.Atom objects (CA atoms)
    """

    try:
        chain = structure[0][chain_id]
    except KeyError:
        raise ValueError(f"Chain {chain_id} not found in {structure.name}")

    ca_atoms = []
    for resnum in residue_numbers:
        for res in chain:
            if res.id[1] == resnum:  # res.id is a tuple: (hetfield, resseq, icode)
                if 'CA' in res:
                    ca_atoms.append(res['CA'])
                else:
                    print(f"Warning: Residue {resnum} does not have a CA atom.")
                break
        else:
            print(f"Warning: Residue {resnum} not found in chain {chain_id}.")

    return ca_atoms

def compute_rmsd(atom_list1, atom_list2):
    """
    Computes RMSD between two lists of atoms.
    """
    coords1 = np.array([atom.get_coord() for atom in atom_list1])
    coords2 = np.array([atom.get_coord() for atom in atom_list2])

    if len(coords1) != len(coords2):
        raise ValueError("Atom lists must be the same length")

    diff = coords1 - coords2
    rmsd = np.sqrt((diff ** 2).sum() / len(coords1))
    return rmsd

def count_initial_gap(seq):
    """
    Counts the number of initial gaps in a sequence.
    """
    count = 0
    for char in seq:
        if char == '-':
            count += 1
        else:
            break
    return count

def align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                            pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2):
    """
    Aligns two PDB files based on MHC and computes RMSD for TCR and peptide chains.
    """
    # Extract chain information
    chain_dict_1 = extract_chain_residue_info(pdbfile_1)
    seq_dict_1 = build_chain_sequences(chain_dict_1)

    chain_dict_2 = extract_chain_residue_info(pdbfile_2)
    seq_dict_2 = build_chain_sequences(chain_dict_2)

    # Align MHC chains
    aligned_seq_mhc_1, aligned_seq_mhc_2, score = align_sequences_middle_gap_penalty(seq_dict_1[mhc_chain_1], seq_dict_2[mhc_chain_2], mode='global')
    # print(f"Aligned MHC Sequences:\n{aligned_seq_mhc_1}\n{aligned_seq_mhc_2}")
    initial_gap_1 = count_initial_gap(aligned_seq_mhc_1)
    initial_gap_2 = count_initial_gap(aligned_seq_mhc_2)
    common_indices = np.array([i for i, (a, b) in enumerate(zip(aligned_seq_mhc_1, aligned_seq_mhc_2)) if a != '-' and b != '-'])
    resnums_mhc_1 = np.array(chain_dict_1[mhc_chain_1]['resnums'])[common_indices - initial_gap_1]
    resnums_mhc_2 = np.array(chain_dict_2[mhc_chain_2]['resnums'])[common_indices - initial_gap_2]
    _, structure1, structure2 = align_chains(pdbfile_1, mhc_chain_1, resnums_mhc_1,
                                                pdbfile_2, mhc_chain_2, resnums_mhc_2)
    mhc_atoms_1 = get_ca_atoms(structure1, mhc_chain_1, resnums_mhc_1)
    mhc_atoms_2 = get_ca_atoms(structure2, mhc_chain_2, resnums_mhc_2)
    print(f'RMSD of MHCs aligning MHCs: {compute_rmsd(mhc_atoms_1, mhc_atoms_2):.2f}')

    # compute RMSD for peptide chain aligned on the MHC
    aligned_seq_pep_1, aligned_seq_pep_2, score = align_sequences_middle_gap_penalty(seq_dict_1[pep_chain_1], seq_dict_2[pep_chain_2], mode='global')
    # print(f"Aligned Peptide Sequences:\n{aligned_seq_pep_1}\n{aligned_seq_pep_2}")
    initial_gap_1 = count_initial_gap(aligned_seq_pep_1)
    initial_gap_2 = count_initial_gap(aligned_seq_pep_2)
    common_indices = np.array([i for i, (a, b) in enumerate(zip(aligned_seq_pep_1, aligned_seq_pep_2)) if a != '-' and b != '-'])
    resnums_pep_1 = np.array(chain_dict_1[pep_chain_1]['resnums'])[common_indices - initial_gap_1]
    resnums_pep_2 = np.array(chain_dict_2[pep_chain_2]['resnums'])[common_indices - initial_gap_2]
    pep_atoms_1 = get_ca_atoms(structure1, pep_chain_1, resnums_pep_1)
    pep_atoms_2 = get_ca_atoms(structure2, pep_chain_2, resnums_pep_2)
    print(f'RMSD of Peptides aligning MHCs: {compute_rmsd(pep_atoms_1, pep_atoms_2):.2f}')
    # print()

    # compute RMSD for TCR alpha and beta chains aligned on the MHC
    aligned_seq_tcra_1, aligned_seq_tcra_2, score = align_sequences_middle_gap_penalty(seq_dict_1[tcra_chain_1], seq_dict_2[tcra_chain_2], mode='global')
    # print(f"Aligned TCR Alpha Sequences:\n{aligned_seq_tcra_1}\n{aligned_seq_tcra_2}")
    initial_gap_1 = count_initial_gap(aligned_seq_tcra_1)
    initial_gap_2 = count_initial_gap(aligned_seq_tcra_2)
    common_indices = np.array([i for i, (a, b) in enumerate(zip(aligned_seq_tcra_1, aligned_seq_tcra_2)) if a != '-' and b != '-'])
    resnums_tcra_1 = np.array(chain_dict_1[tcra_chain_1]['resnums'])[common_indices - initial_gap_1]
    resnums_tcra_2 = np.array(chain_dict_2[tcra_chain_2]['resnums'])[common_indices - initial_gap_2]

    aligned_seq_tcrb_1, aligned_seq_tcrb_2, score = align_sequences_middle_gap_penalty(seq_dict_1[tcrb_chain_1], seq_dict_2[tcrb_chain_2], mode='global')
    # print(f"Aligned TCR Beta Sequences:\n{aligned_seq_tcrb_1}\n{aligned_seq_tcrb_2}")
    initial_gap_1 = count_initial_gap(aligned_seq_tcrb_1)
    initial_gap_2 = count_initial_gap(aligned_seq_tcrb_2)
    common_indices = np.array([i for i, (a, b) in enumerate(zip(aligned_seq_tcrb_1, aligned_seq_tcrb_2)) if a != '-' and b != '-'])
    resnums_tcrb_1 = np.array(chain_dict_1[tcrb_chain_1]['resnums'])[common_indices - initial_gap_1]
    resnums_tcrb_2 = np.array(chain_dict_2[tcrb_chain_2]['resnums'])[common_indices - initial_gap_2]
    tcr_atoms_1 = get_ca_atoms(structure1, tcra_chain_1, resnums_tcra_1) + get_ca_atoms(structure1, tcrb_chain_1, resnums_tcrb_1)
    tcr_atoms_2 = get_ca_atoms(structure2, tcra_chain_2, resnums_tcra_2) + get_ca_atoms(structure2, tcrb_chain_2, resnums_tcrb_2)
    print(f'RMSD of TCRs aligning MHCs: {compute_rmsd(tcr_atoms_1, tcr_atoms_2):.2f}')
    # print()

    # also do a full alignment of the whole system
    super_imposer, structure2 = align_by_atoms(structure2, mhc_atoms_1 + pep_atoms_1 + tcr_atoms_1, mhc_atoms_2 + pep_atoms_2 + tcr_atoms_2)
    print(f"RMSD of full alignment: {super_imposer.rms:.3f}")
    print()
    print()


if __name__ == '__main__':

    print()
    print('tcr1 - no template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/mskcc/pdbs/5d2n-filtered.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'I'
    tcra_chain_1 = 'D'
    tcrb_chain_1 = 'E'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/mskcc/pdbs/fold_tcr1_no_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)

    print('tcr1 - yes template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/mskcc/pdbs/5d2n-filtered.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'I'
    tcra_chain_1 = 'D'
    tcrb_chain_1 = 'E'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/mskcc/pdbs/fold_tcr1_yes_template_model_1.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)

    print('nyeso - no template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/nyeso/pdbs/2bnr.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'C'
    tcra_chain_1 = 'D'
    tcrb_chain_1 = 'E'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/nyeso/pdbs/fold_nyeso_no_template_model_1.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('nyeso - yes template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/nyeso/pdbs/2bnr.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'C'
    tcra_chain_1 = 'D'
    tcrb_chain_1 = 'E'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/nyeso/pdbs/fold_nyeso_yes_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('nyeso 9v - no template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/nyeso/pdbs/2bnq.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'C'
    tcra_chain_1 = 'D'
    tcrb_chain_1 = 'E'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/nyeso/pdbs/fold_nyeso_9v_no_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('nyeso 9v - yes template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/nyeso/pdbs/2bnq.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'C'
    tcra_chain_1 = 'D'
    tcrb_chain_1 = 'E'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/nyeso/pdbs/fold_nyeso_9v_yes_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('tax - no template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/tax/pdbs/1ao7.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'C'
    tcra_chain_1 = 'D'
    tcrb_chain_1 = 'E'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/tax/pdbs/fold_tax_no_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)

    print('tax - yes template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/tax/pdbs/1ao7.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'C'
    tcra_chain_1 = 'D'
    tcrb_chain_1 = 'E'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/tax/pdbs/fold_tax_yes_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)

    print('tax 8a - no template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/tax/pdbs/1qsf.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'C'
    tcra_chain_1 = 'D'
    tcrb_chain_1 = 'E'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/tax/pdbs/fold_tax_8a_no_template_model_1.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)

    print('tax 8a - yes template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/tax/pdbs/1qsf.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'C'
    tcra_chain_1 = 'D'
    tcrb_chain_1 = 'E'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/tax/pdbs/fold_tax_8a_yes_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)

    print('hsiue - no template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/hsiue_et_al/pdbs/6w51-filtered.pdb'
    mhc_chain_1 = 'D'
    pep_chain_1 = 'F'
    tcra_chain_1 = 'O'
    tcrb_chain_1 = 'P'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/hsiue_et_al/pdbs/fold_hsiue_no_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'C'
    tcra_chain_2 = 'D'
    tcrb_chain_2 = 'E'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('hsiue - yes template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/hsiue_et_al/pdbs/6w51-filtered.pdb'
    mhc_chain_1 = 'D'
    pep_chain_1 = 'F'
    tcra_chain_1 = 'O'
    tcrb_chain_1 = 'P'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/hsiue_et_al/pdbs/fold_hsiue_yes_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'C'
    tcra_chain_2 = 'D'
    tcrb_chain_2 = 'E'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('mage - no template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/pdbs/5brz.pdb.human.MH1.A-01.A.C.DE.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'B'
    tcra_chain_1 = 'C'
    tcrb_chain_1 = 'D'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/pdbs/fold_mage_no_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('mage - yes template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/pdbs/5brz.pdb.human.MH1.A-01.A.C.DE.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'B'
    tcra_chain_1 = 'C'
    tcrb_chain_1 = 'D'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/pdbs/fold_mage_yes_template_seed_2_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('titin - no template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/pdbs/5bs0.pdb.human.MH1.A-01.A.C.DE.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'B'
    tcra_chain_1 = 'C'
    tcrb_chain_1 = 'D'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/pdbs/fold_titin_no_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('titin - yes template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/pdbs/5bs0.pdb.human.MH1.A-01.A.C.DE.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'B'
    tcra_chain_1 = 'C'
    tcrb_chain_1 = 'D'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/pdbs/fold_titin_yes_template_model_1.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('ebv - no template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/pdbs/3mv7.pdb.human.MH1.B-35.A.C.DE.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'B'
    tcra_chain_1 = 'C'
    tcrb_chain_1 = 'D'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/pdbs/fold_ebv_no_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)

    print('ebv - yes template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/pdbs/3mv7.pdb.human.MH1.B-35.A.C.DE.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'B'
    tcra_chain_1 = 'C'
    tcrb_chain_1 = 'D'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/pdbs/fold_ebv_yes_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)


    print('ebv q5 - no template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/pdbs/4prp.pdb.human.MH1.B-35.A.C.DE.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'B'
    tcra_chain_1 = 'C'
    tcrb_chain_1 = 'D'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/pdbs/fold_ebv_q5_no_template_model_0.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)

    print('ebv q5 - yes template')

    pdbfile_1 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/pdbs/4prp.pdb.human.MH1.B-35.A.C.DE.pdb'
    mhc_chain_1 = 'A'
    pep_chain_1 = 'B'
    tcra_chain_1 = 'C'
    tcrb_chain_1 = 'D'

    pdbfile_2 = '/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/pdbs/fold_ebv_q5_yes_template_model_1.pdb'
    mhc_chain_2 = 'A'
    pep_chain_2 = 'B'
    tcra_chain_2 = 'C'
    tcrb_chain_2 = 'D'

    align_on_mhc_and_compute_rmsds(pdbfile_1, mhc_chain_1, pep_chain_1, tcra_chain_1, tcrb_chain_1,
                                pdbfile_2, mhc_chain_2, pep_chain_2, tcra_chain_2, tcrb_chain_2)

