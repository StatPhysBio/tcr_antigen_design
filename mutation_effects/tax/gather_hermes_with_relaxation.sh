
model_version_list='hermes_py_000 hermes_py_050 hermes_py_000_ft_skempi_no_tcrpmhc_ddg_bi hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi'

pdbs='1ao7 1qse 1qsf'

metrics='pnE pnlogp'

for model_version in $model_version_list
    do

    for metric in $metrics
        do

        # for pdb in $pdbs
        #     do

        #     python ../../src/gather_score_with_relaxation_results.py \
        #                     --experiment_dir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/tax \
        #                     --csv_filename 'tax_peptide_kd_'$pdb'.csv' \
        #                     --pdbid $pdb \
        #                     --model_version $model_version
            
        #     python -u pretty_plots.py \
        #                 --system tax \
        #                 --system_name_in_csv_file 'tax_peptide_kd_'$pdb'_with_relaxation' \
        #                 --target_column'=-log10(Kd)' \
        #                 --prediction_column $metric \
        #                 --model hermes \
        #                 --model_instance $model_version \
        #                 --use_mt_structure 0 \
        #                 --show_wt_lines both_from_df
            
        # done

        # # all pdbs
        # python ../../src/gather_score_with_relaxation_results.py \
        #                 --experiment_dir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/tax \
        #                 --csv_filename 'tax_peptide_kd_averaged.csv' \
        #                 --pdbid all \
        #                 --model_version $model_version
        
        # python -u pretty_plots.py \
        #             --system tax \
        #             --system_name_in_csv_file 'tax_peptide_kd_averaged_with_relaxation' \
        #             --target_column'=-log10(Kd)' \
        #             --prediction_column $metric \
        #             --model hermes \
        #             --model_instance $model_version \
        #             --use_mt_structure 0 \
        #             --show_wt_lines both_from_df

        # # closest pdbs
        # python ../../src/gather_score_with_relaxation_results.py \
        #                 --experiment_dir ./ \
        #                 --csv_filename 'tax_peptide_kd_closest.csv' \
        #                 --pdbid from_csv \
        #                 --model_version $model_version \
        #                 --pdb_column pdb
        
        # python -u pretty_plots.py \
        #             --base_dir ../ \
        #             --system tax \
        #             --system_name_in_csv_file 'tax_peptide_kd_closest_with_relaxation' \
        #             --target_column'=-log10(Kd)' \
        #             --prediction_column $metric \
        #             --model hermes \
        #             --model_instance $model_version \
        #             --use_mt_structure 0 \
        #             --show_wt_lines both_from_df


        python ../src/gather_score_with_relaxation_results.py \
                        --experiment_dir ./ \
                        --csv_filename 'tax_peptide_kd_closest.csv' \
                        --pdbid from_csv \
                        --model_version $model_version \
                        --pdb_column pdb \
                        --use_min_rosetta_energy_instead_of_full_average 1
        
        python -u pretty_plots.py \
                    --base_dir ../ \
                    --system tax \
                    --system_name_in_csv_file 'tax_peptide_kd_closest_with_relaxation_min_energy' \
                    --target_column'=-log10(Kd)' \
                    --prediction_column $metric \
                    --model hermes \
                    --model_instance $model_version \
                    --use_mt_structure 0 \
                    --show_wt_lines both_from_df
        

        python ../src/gather_score_with_relaxation_results.py \
                        --experiment_dir ./ \
                        --csv_filename 'tax_peptide_kd_closest.csv' \
                        --pdbid from_csv \
                        --model_version $model_version \
                        --pdb_column pdb \
                        --use_min_rosetta_energy_runs_but_compute_mean 1
        
        python -u pretty_plots.py \
                    --base_dir ../ \
                    --system tax \
                    --system_name_in_csv_file 'tax_peptide_kd_closest_with_relaxation_mean_but_min_energy_runs' \
                    --target_column'=-log10(Kd)' \
                    --prediction_column $metric \
                    --model hermes \
                    --model_instance $model_version \
                    --use_mt_structure 0 \
                    --show_wt_lines both_from_df
    
    done

done
