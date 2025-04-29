

python -u pretty_plots.py \
            --base_dir ../ \
            --system nyeso \
            --system_name_in_csv_file nyeso_peptide_kd_closest \
            --target_column'=-log10(Kd)' \
            --model tapir \
            --model_instance tapir \
            --show_wt_lines both_from_df


