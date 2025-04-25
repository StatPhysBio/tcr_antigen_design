


base_dir='./'
output_dir=$base_dir'results/'

sub_matrices='blosum62 luksza_cross_reactivity luksza_cross_reactivity_without_d'

for sub_matrix in $sub_matrices
    do

    echo $sub_matrix

    python -u ../src/score_with_sub_matrix.py \
                --input_csv_filepath 'hsiue_et_al_H2_sat_mut.csv' \
                --output_csv_filepath $output_dir'/'$sub_matrix'/hsiue_et_al_H2_sat_mut-'$sub_matrix'-use_mt_structure=0.csv' \
                --wt_peptide_column sequence \
                --mt_peptide_column wt_seq \
                --substitution_matrix $sub_matrix

    python -u pretty_plots.py \
                --system hsiue_et_al \
                --system_name_in_csv_file hsiue_et_al_H2_sat_mut \
                --target_column'=IFN_gamma (pg/ml)' \
                --model $sub_matrix \
                --model_instance $sub_matrix \
                --show_wt_lines both_from_df

done