


base_dir='./'
output_dir=$base_dir'results/'


python -u ../../src/score_with_sub_matrix.py \
            --input_csv_filepath 'hsiue_et_al_H2_sat_mut.csv' \
            --output_csv_filepath $output_dir'/blosum62/hsiue_et_al_H2_sat_mut-blosum62-use_mt_structure=0.csv' \
            --wt_peptide_column sequence \
            --mt_peptide_column wt_seq

python -u pretty_plots.py \
            --system hsiue_et_al \
            --system_name_in_csv_file hsiue_et_al_H2_sat_mut \
            --target_column'=IFN_gamma (pg/ml)' \
            --model blosum62 \
            --model_instance blosum62 \
            --show_wt_lines both_from_df

