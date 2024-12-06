
tcrs='1 2 3 4 5 6 7'

use_mt_structure='0'

model_version_list='tcrdock tcrdock_no_nearby_templates'

for tcr in $tcrs
    do
    for model_version in $model_version_list
        do

        python -u pretty_plots.py \
                        --system mskcc \
                        --system_name_in_csv_file 'mskcc_tcr'$tcr'_ec50_sat_mut' \
                        --target_column'=- delta log_ec50_M' \
                        --model $model_version \
                        --model_instance $model_version \
                        --num_seq_per_target 10 \
                        --use_mt_structure $use_mt_structure \
                        --show_wt_lines target \
                        --wt_value_pred 0.0
                        

    done
done

