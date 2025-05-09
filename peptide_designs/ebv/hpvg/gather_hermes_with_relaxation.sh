


python ../../src/gather_score_with_relaxation_results.py \
                --csv_filename hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv \
                --pdbid from_csv \
                --model_version hermes_py_000 \
                --pdb_column pdb_for_hermes_scoring \
                --seq_column peptide \
                --output_in_inputfile 1 \
                --use_max_instead_of_mean 1 # more closely approximates the real pnE that was computed during annealing

# python ../../src/gather_score_with_relaxation_results.py \
#                 --csv_filename hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv \
#                 --pdbid from_csv \
#                 --model_version hermes_py_050 \
#                 --pdb_column pdb_for_hermes_scoring \
#                 --seq_column peptide \
#                 --output_in_inputfile 1
        
