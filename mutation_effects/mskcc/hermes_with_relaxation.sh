

# python ../../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py \
#                             --model_version hermes_py_000 \
#                             --pdb TCR2_T00000_A0201_NLVPMVATV_0_model_1_model_2_ptm_ft4 \
#                             --pdbdir ./pdbs \
#                             --chain A \
#                             --peptide_resnum_start 376 \
#                             --sequence NLVPMVATVE \
#                             --num_repeats 1 \
#                             --output_dir ./test_output \
#                             --verbose 0 \
#                             --job 0

python ../../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py \
                            --model_version hermes_py_000 \
                            --pdb 5d2n-filtered \
                            --pdbdir ./pdbs \
                            --chain I \
                            --peptide_resnum_start 1 \
                            --sequence NLVPMVATVE \
                            --num_repeats 1 \
                            --output_dir ./test_output \
                            --verbose 0 \
                            --job 0
