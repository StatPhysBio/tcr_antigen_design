
model_version_list='hermes_py_000 hermes_py_050'

pdbid_list=(
    5d2n-filtered
    af3_tcr2
    af3_tcr3
    af3_tcr4
    af3_tcr5
    af3_tcr6
    af3_tcr7
)

metrics='pnE'

for model_version in $model_version_list
    do

    for tcr in $(seq 1 7)
        do

        pdbid=${pdbid_list[$tcr-1]}

        for metric in $metrics
            do

            echo $tcr $pdbid $metric $model_version

            python ../src/gather_score_with_relaxation_results.py \
                            --experiment_dir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc \
                            --csv_filename 'mskcc_tcr'$tcr'_ec50_sat_mut_af3_no_template.csv' \
                            --pdbid from_csv \
                            --pdb_column wt_pdb \
                            --model_version $model_version
            
            python -u pretty_plots.py \
                        --system mskcc \
                        --system_name_in_csv_file 'mskcc_tcr'$tcr'_ec50_sat_mut_af3_no_template_with_relaxation' \
                        --target_column'=-log10(EC50)' \
                        --prediction_column $metric \
                        --model hermes \
                        --model_instance $model_version \
                        --use_mt_structure 0 \
                        --show_wt_lines both_from_df

        
            # python ../src/gather_score_with_relaxation_results.py \
            #                 --experiment_dir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc \
            #                 --csv_filename 'mskcc_tcr'$tcr'_ec50_sat_mut_af3.csv' \
            #                 --pdbid $pdbid \
            #                 --model_version $model_version \
            #                 --use_min_rosetta_energy_instead_of_full_average 1
            
            # python -u pretty_plots.py \
            #             --system mskcc \
            #             --system_name_in_csv_file 'mskcc_tcr'$tcr'_ec50_sat_mut_af3_with_relaxation_min_energy' \
            #             --target_column'=-log10(EC50)' \
            #             --prediction_column $metric \
            #             --model hermes \
            #             --model_instance $model_version \
            #             --use_mt_structure 0 \
            #             --show_wt_lines both_from_df

            # python ../src/gather_score_with_relaxation_results.py \
            #                 --experiment_dir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc \
            #                 --csv_filename 'mskcc_tcr'$tcr'_ec50_sat_mut_af3.csv' \
            #                 --pdbid $pdbid \
            #                 --model_version $model_version \
            #                 --use_min_rosetta_energy_runs_but_compute_mean 1
            
            # python -u pretty_plots.py \
            #             --system mskcc \
            #             --system_name_in_csv_file 'mskcc_tcr'$tcr'_ec50_sat_mut_af3_with_relaxation_mean_but_min_energy_runs' \
            #             --target_column'=-log10(EC50)' \
            #             --prediction_column $metric \
            #             --model hermes \
            #             --model_instance $model_version \
            #             --use_mt_structure 0 \
            #             --show_wt_lines both_from_df

        done
    
    done

    for metric in $metrics
        do

        python ../src/gather_score_with_relaxation_results.py \
                        --experiment_dir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc \
                        --csv_filename mskcc_tcr1_ec50_sat_mut_af3_yes_template.csv \
                        --pdbid from_csv \
                        --pdb_column wt_pdb \
                        --model_version $model_version
        
        python -u pretty_plots.py \
                    --system mskcc \
                    --system_name_in_csv_file mskcc_tcr1_ec50_sat_mut_af3_yes_template_with_relaxation \
                    --target_column'=-log10(EC50)' \
                    --prediction_column $metric \
                    --model hermes \
                    --model_instance $model_version \
                    --use_mt_structure 0 \
                    --show_wt_lines both_from_df
        
        done
done
