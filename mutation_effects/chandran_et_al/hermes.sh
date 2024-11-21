

hermes_dir='/gscratch/spe/gvisan01/hermes/'

tcrs='3 4'

model_version_list='hermes_py_000 hermes_py_050 hermes_py_000_ft_skempi_no_tcrpmhc_ddg_bi hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi'

use_mt_structure='0'

base_dir='./'
output_dir=$base_dir'results/'

for tcr in $tcrs
    do

    echo $tcr

    for model_version in $model_version_list
        do

        echo $model_version

        # python -u $hermes_dir'mutation_effect_prediction_with_hermes.py' \
        #                     --model_version $model_version \
        #                     --csv_file $base_dir'chandran_et_al_peptide_A_and_G_scans_TCR'$tcr'.csv' \
        #                     --folder_with_pdbs ./pdbs \
        #                     --output_dir $output_dir \
        #                     --wt_pdb_column wt_pdb \
        #                     --mutant_column mutant \
        #                     --mutant_chain_column mutant_chain \
        #                     --use_mt_structure $use_mt_structure

        python -u pretty_plots.py \
                    --system chandran_et_al \
                    --system_name_in_csv_file 'chandran_et_al_peptide_A_and_G_scans_TCR'$tcr \
                    --target_column'=TFN_alpha_mt/TFN_alpha_wt' \
                    --model hermes \
                    --model_instance $model_version \
                    --use_mt_structure $use_mt_structure \
                    --show_wt_lines both \
                    --wt_value_target 100.0
    
    done
done

