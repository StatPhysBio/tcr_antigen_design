


python -u ../../../mutation_effects/src/get_score_of_requested_peptides_fixed_structure.py \
            --model_version hermes_py_000 \
            --pdbdir ./pdbs \
            --output_csv_filepath hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv \
            --csv_filepath hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv \
            --pdb_column pdb_for_hermes_scoring \
            --chain_column chain_for_hermes_scoring \
            --peptide_column peptide \
            --peptide_resnum_start 181

python -u ../../../mutation_effects/src/get_score_of_requested_peptides_fixed_structure.py \
            --model_version hermes_py_050 \
            --pdbdir ./pdbs \
            --output_csv_filepath hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv \
            --csv_filepath hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv \
            --pdb_column pdb_for_hermes_scoring \
            --chain_column chain_for_hermes_scoring \
            --peptide_column peptide \
            --peptide_resnum_start 181



