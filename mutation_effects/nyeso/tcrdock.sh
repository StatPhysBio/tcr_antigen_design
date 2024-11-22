


pdbs='2bnr 2bnq'


use_mt_structure='0'

model_version_list='tcrdock tcrdock_no_nearby_templates'


for model_version in $model_version_list
    do

    python -u pretty_plots.py \
                    --system nyeso \
                    --system_name_in_csv_file nyeso_peptide_kd \
                    --target_column'=-log10(Kd)' \
                    --model $model_version \
                    --model_instance $model_version \
                    --num_seq_per_target 10 \
                    --use_mt_structure $use_mt_structure \
                    

done


