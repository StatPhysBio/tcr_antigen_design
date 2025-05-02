
import gemmi

def remove_symmetry_and_save_asu(input_pdb, output_pdb):
    # Read structure
    st = gemmi.read_structure(input_pdb)
    
    # Only keep the first model
    model = st[0]

    # Remove symmetry-expanded chains (i.e., keep only chains with identity operator)
    # Often, only chains with identity operator belong to the asymmetric unit
    # However, gemmi read_structure on a PDB (not mmCIF) already gives you just ASU

    # Write out the cleaned structure
    st.write_pdb(output_pdb)
    print(f"Saved asymmetric unit to {output_pdb}")

def build_biological_assembly(input_pdb, output_pdb, assembly_id='1'):
    doc = gemmi.cif.read_file(input_pdb)  # also works with mmCIF
    st = gemmi.make_assembly_from_mmcif(doc, assembly_id)
    st.write_pdb(output_pdb)
    print(f"Wrote biological assembly {assembly_id} to {output_pdb}")

build_biological_assembly('5d2n.pdb', '5d2n_asu.pdb')
