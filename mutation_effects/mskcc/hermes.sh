

hermes_dir='/gscratch/spe/gvisan01/hermes/'

tcrs='2 3 4 5 6 7' # 1 2 3 4 5 6 7

model_version_list='hermes_py_000 hermes_py_050'

use_mt_structure='0'

base_dir='./'
output_dir=$base_dir'results/'

for tcr in $tcrs
    do

    echo $tcr

    for model_version in $model_version_list
        do

        echo $model_version

        python -u $hermes_dir'mutation_effect_prediction_with_hermes.py' \
                            --model_version $model_version \
                            --csv_file $base_dir'mskcc_tcr'$tcr'_ec50_sat_mut_af3.csv' \
                            --folder_with_pdbs ./pdbs \
                            --output_dir $output_dir \
                            --wt_pdb_column wt_pdb \
                            --mutant_column mutant \
                            --mutant_chain_column mutant_chain \
                            --use_mt_structure $use_mt_structure

        python -u pretty_plots.py \
                    --system mskcc \
                    --system_name_in_csv_file 'mskcc_tcr'$tcr'_ec50_sat_mut_af3' \
                    --target_column'=- delta log_ec50_M' \
                    --model hermes \
                    --model_instance $model_version \
                    --use_mt_structure $use_mt_structure \
                    --show_wt_lines both
    
    done
done

# '=- delta log_ec50_M' \
# '=<Theta21>' \
