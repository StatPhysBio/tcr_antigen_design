
scoring_functions_dir='/gscratch/spe/gvisan01/hermes/baselines/'

tcrs='8'

matrix_list='blosum62 neg_abs_diff_vdw_radius'
subtract_wt_value_list='1'
combine_fn_list='sum'

base_dir='./'
output_dir=$base_dir'results/'

for tcr in $tcrs
    do

    echo $tcr

    for matrix in $matrix_list
        do
        for subtract_wt_value in $subtract_wt_value_list
            do
            for combine_fn in $combine_fn_list
                do

                echo $matrix $subtract_wt_value $combine_fn

                matrix_plus_stuff=$matrix'__'$subtract_wt_value'__'$combine_fn

                python -u $scoring_functions_dir'score_with_substitution_matrix.py' \
                            --csv_file $base_dir'hhatv4_tcr'$tcr'_ec50_sat_mut.csv' \
                            --mutant_column mutant \
                            --mutant_split_symbol"=|" \
                            --output_dir $output_dir \
                            --substitution_matrix $matrix \
                            --subtract_wt_value $subtract_wt_value \
                            --combine_multiple_mutations_with $combine_fn \

                python -u pretty_plots.py \
                            --system hhatv4 \
                            --system_name_in_csv_file 'hhatv4_tcr'$tcr'_ec50_sat_mut' \
                            --target_column'=- delta log_ec50_M' \
                            --model $matrix \
                            --model_instance $matrix_plus_stuff \
                            --show_wt_lines both

            done
        done
    done
done