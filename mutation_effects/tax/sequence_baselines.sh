

python -u pretty_plots.py \
            --base_dir ../ \
            --system tax \
            --system_name_in_csv_file tax_peptide_kd_closest \
            --target_column'=-log10(Kd)' \
            --model tapir \
            --model_instance tapir \
            --show_wt_lines both_from_df


