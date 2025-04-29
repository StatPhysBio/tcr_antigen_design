
python -u ../src/score_peptides_with_esmif.py \
            --pdbdir ./pdbs \
            --output_csv_filepath ./esmif/nyeso_peptide_kd_closest-esmif-use_mt_structure=0.csv \
            --csv_filepath nyeso_peptide_kd_closest.csv \
            --pdb_column pdb \
            --chain_column chain \
            --peptide_column sequence
