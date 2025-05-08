
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

    from Bio.PDB import PDBParser, is_aa

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
        prev_num = resnums[0] - 1
        
        for resname, resnum in zip(resnames, resnums):
            gap = resnum - prev_num - 1  # How many residues missing?
            if gap > 0:
                sequence += '-' * gap  # Insert gaps
            sequence += resname
            prev_num = resnum
        
        seq_dict[chain_id] = sequence
    
    return seq_dict


seq_dict = build_chain_sequences(extract_chain_residue_info('6w51-filtered.pdb'))
print(seq_dict)

