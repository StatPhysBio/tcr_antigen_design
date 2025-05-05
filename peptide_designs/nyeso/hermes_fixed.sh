
model_version_list='hermes_py_000 hermes_py_050'

for model_version in $model_version_list
    do

    python -u ../../mutation_effects/src/get_score_of_requested_peptides_fixed_structure.py \
                --model_version $model_version \
                --pdbdir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/pdbs \
                --output_csv_filepath 'hermes_scores/'$model_version'/selected_designs_for_hermes_scoring__'$model_version'__fixed' \
                --csv_filepath selected_designs_for_hermes_scoring.csv \
                --pdb_column pdb \
                --chain_column chain \
                --peptide_column sequence \
                --peptide_resnum_start 1

done