
model_version_list='hermes_py_000 hermes_py_050'

pdbid_list=(
    5d2n-filtered
    TCR2_T00000_A0201_NLVPMVATV_0_model_1_model_2_ptm_ft4
    TCR3_T00000_A0201_NLVPMVATV_0_model_1_model_2_ptm_ft4
    TCR4_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4
    TCR5_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4
    TCR6_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4
    TCR7_T00000_B2705_GRLKALCQR_0_model_1_model_2_ptm_ft4
)

metrics='pnE pnlogp'

for model_version in $model_version_list
    do

    for tcr in $(seq 1 1)
        do

        pdbid=${pdbid_list[$tcr-1]}

        for metric in $metrics
            do


            python ../../src/gather_score_with_relaxation_results.py \
                            --experiment_dir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc \
                            --csv_filename 'mskcc_tcr'$tcr'_ec50_sat_mut.csv' \
                            --pdbid $pdbid \
                            --model_version $model_version
            
            python -u pretty_plots.py \
                        --system mskcc \
                        --system_name_in_csv_file 'mskcc_tcr'$tcr'_ec50_sat_mut_with_relaxation' \
                        --target_column'=- delta log_ec50_M' \
                        --prediction_column $metric \
                        --model hermes \
                        --model_instance $model_version \
                        --use_mt_structure 0 \
                        --show_wt_lines both_from_df


        done
    
    done

done
