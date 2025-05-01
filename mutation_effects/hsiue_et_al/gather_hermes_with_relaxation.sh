
model_version_list='hermes_py_000 hermes_py_050'

metrics='pnE pnlogp'

for model_version in $model_version_list
    do

    for metric in $metrics
        do


        python ../src/gather_score_with_relaxation_results.py \
                        --experiment_dir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/hsiue_et_al \
                        --csv_filename 'hsiue_et_al_H2_sat_mut.csv' \
                        --pdbid 6w51-filtered \
                        --model_version $model_version \
                        --use_min_rosetta_energy_runs_but_compute_mean 1 \
                        # --use_min_rosetta_energy_instead_of_full_average 1
        
        python -u pretty_plots.py \
                    --system hsiue_et_al \
                    --system_name_in_csv_file 'hsiue_et_al_H2_sat_mut_with_relaxation_mean_but_min_energy_runs' \
                    --target_column'=IFN_gamma (pg/ml)' \
                    --prediction_column $metric \
                    --model hermes \
                    --model_instance $model_version \
                    --use_mt_structure 0 \
                    --show_wt_lines both_from_df


    done

done
