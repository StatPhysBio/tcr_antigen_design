
base_dir='/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/'

python -u ../src/score_peptides_with_esmif.py \
            --pdbdir $base_dir'pdbs' \
            --output_csv_filepath $base_dir'results/esmif/nyeso_peptide_kd_closest-esmif-use_mt_structure=0.csv' \
            --csv_filepath $base_dir'nyeso_peptide_kd_closest.csv' \
            --pdb_column pdb \
            --chain_column chain \
            --peptide_column sequence
