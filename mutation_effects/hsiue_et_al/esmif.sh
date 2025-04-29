
python -u ../src/score_peptides_with_esmif.py \
            --pdbdir ./pdbs \
            --output_csv_filepath ./esmif/hsiue_et_al_H2_sat_mut-esmif-use_mt_structure=0.csv \
            --csv_filepath hsiue_et_al_H2_sat_mut.csv \
            --pdb_column pdb \
            --chain_column chain \
            --peptide_column sequence
