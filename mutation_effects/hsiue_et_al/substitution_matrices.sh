
scoring_functions_dir='/gscratch/spe/gvisan01/hermes/baselines/'

matrix_list='blosum62 neg_abs_diff_vdw_radius'
subtract_wt_value_list='1'
combine_fn_list='sum'

base_dir='./'
output_dir=$base_dir'results/'

for matrix in $matrix_list
    do
    for subtract_wt_value in $subtract_wt_value_list
        do
        for combine_fn in $combine_fn_list
            do

            echo $matrix $subtract_wt_value $combine_fn

            matrix_plus_stuff=$matrix'__'$subtract_wt_value'__'$combine_fn

            python -u $scoring_functions_dir'score_with_substitution_matrix.py' \
                        --csv_file $base_dir'hsiue_et_al_H2_sat_mut.csv' \
                        --mutant_column mutant \
                        --mutant_split_symbol"=|" \
                        --output_dir $output_dir \
                        --substitution_matrix $matrix \
                        --subtract_wt_value $subtract_wt_value \
                        --combine_multiple_mutations_with $combine_fn \

            python -u pretty_plots.py \
                        --system hsiue_et_al \
                        --system_name_in_csv_file 'hsiue_et_al_H2_sat_mut' \
                        --target_column'=IFN_gamma (pg/ml)' \
                        --model $matrix \
                        --model_instance $matrix_plus_stuff \
                        --show_wt_lines both \
                        --wt_value_target 10596.0 \
                        --wt_value_pred 0.0

        done
    done
done
