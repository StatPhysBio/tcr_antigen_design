


base_dir='./'
output_dir=$base_dir'results/'

tcrs='1 2 3 4 5 6 7'

for tcr in $tcrs
    do

    echo $tcr

    python -u ../../src/score_with_sub_matrix.py \
                --input_csv_filepath 'mskcc_tcr'$tcr'_ec50_sat_mut_af3.csv' \
                --output_csv_filepath $output_dir'/blosum62/mskcc_tcr'$tcr'_ec50_sat_mut_af3-blosum62-use_mt_structure=0.csv' \
                --wt_peptide_column sequence \
                --mt_peptide_column wt_seq

    python -u pretty_plots.py \
                --system mskcc \
                --system_name_in_csv_file 'mskcc_tcr'$tcr'_ec50_sat_mut_af3' \
                --target_column'=- delta log_ec50_M' \
                --model blosum62 \
                --model_instance blosum62 \
                --show_wt_lines both_from_df

done
