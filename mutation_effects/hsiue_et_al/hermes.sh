

hermes_dir='/gscratch/spe/gvisan01/hermes/'

model_version_list='hermes_py_000 hermes_py_050 hermes_py_000_ft_skempi_ddg_bi hermes_py_050_ft_skempi_ddg_bi'

use_mt_structure='0'

base_dir='./'
output_dir=$base_dir'results/'

for model_version in $model_version_list
    do

    echo $model_version

    python -u $hermes_dir'mutation_effect_prediction_with_hermes.py' \
                        --model_version $model_version \
                        --csv_file $base_dir'hsiue_et_al_H2_sat_mut.csv' \
                        --folder_with_pdbs ./pdbs \
                        --output_dir $output_dir \
                        --wt_pdb_column wt_pdb \
                        --mutant_column mutant \
                        --mutant_chain_column mutant_chain \
                        --use_mt_structure $use_mt_structure

    python -u pretty_plots.py \
                --system hsiue_et_al \
                --system_name_in_csv_file 'hsiue_et_al_H2_sat_mut' \
                --target_column'=IFN_gamma (pg/ml)' \
                --model hermes \
                --model_instance $model_version \
                --use_mt_structure $use_mt_structure \
                --show_wt_lines both \
                --wt_value_target 10596.0 \
                --wt_value_pred 0.0

done
