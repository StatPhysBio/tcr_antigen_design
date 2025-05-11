
model_version_list='hermes_py_000 hermes_py_050'

for model_version in $model_version_list
    do

    python ../src/gather_score_with_relaxation_results.py \
                    --csv_filename selected_designs_for_hermes_scoring.csv \
                    --pdbid from_csv \
                    --model_version $model_version \
                    --pdb_column pdb \
                    --use_max_instead_of_mean 0 # if 1, more closely approximates the real pnE that was computed during annealing (maybe)
        
done
