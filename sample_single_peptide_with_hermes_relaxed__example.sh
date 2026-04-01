

python sample_single_peptide_with_hermes_relaxed.py \
    --hermes_dir /gscratch/stf/gvisan01/hermes/ \
    -m hermes_py_000 \
    -p /gscratch/stf/gvisan01/tcr_antigen_design/peptide_designs/magea3_and_titin/magea3/pdbs/5brz.pdb.human.MH1.A-01.A.C.DE.pdb \
    -o test_relax_0.txt \
    -c B \
    --initial_sequence AAAAAAAAA \
    --schedule_name multiplicative \
    --iters 3 \
    --save_metric_history \
    --write_pdb
