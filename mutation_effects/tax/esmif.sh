
python -u ../src/score_peptides_with_esmif.py \
            --pdbdir ./pdbs \
            --output_csv_filepath ./esmif/tax_peptide_kd_closest-esmif-use_mt_structure=0.csv \
            --csv_filepath tax_peptide_kd_closest.csv \
            --pdb_column pdb \
            --chain_column chain \
            --peptide_column sequence
