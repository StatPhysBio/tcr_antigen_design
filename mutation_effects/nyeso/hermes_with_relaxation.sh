
python ../../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py \
            --model_version hermes_py_000 \
            --pdb 2bnq \
            --pdbdir ./pdbs \
            --chain C \
            --sequence SLLMYITQV \
            --num_repeats 3 \
            --output_dir ./results/hermes_py_000/with_relaxation/ \
            --verbose 0

