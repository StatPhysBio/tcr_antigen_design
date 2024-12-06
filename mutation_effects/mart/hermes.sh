

hermes_dir='/gscratch/spe/gvisan01/hermes/'

pdbs='3qdg 6am5 6amu'

model_version_list='hermes_py_000 hermes_py_050 hermes_py_000_ft_skempi_no_tcrpmhc_ddg_bi hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi'

use_mt_structure='0'

base_dir='./'
output_dir=$base_dir'results/'

for pdb in $pdbs
    do

    echo $pdb

    for model_version in $model_version_list
        do

        echo $model_version

        # python -u $hermes_dir'mutation_effect_prediction_with_hermes.py' \
        #                     --model_version $model_version \
        #                     --csv_file $base_dir'mart_peptide_kd_'$pdb'.csv' \
        #                     --folder_with_pdbs ./pdbs \
        #                     --output_dir $output_dir \
        #                     --wt_pdb_column pdb \
        #                     --mutant_column mutants \
        #                     --mutant_chain_column chain \
        #                     --mutant_split_symbol '|' \
        #                     --use_mt_structure $use_mt_structure

        python -u pretty_plots.py \
                    --system mart \
                    --system_name_in_csv_file 'mart_peptide_kd_'$pdb \
                    --target_column'=neg_delta_log10_Kd' \
                    --model hermes \
                    --model_instance $model_version \
                    --use_mt_structure $use_mt_structure \
                    --show_wt_lines both \
                    --wt_value_target 0.0 \
                    --wt_value_pred 0.0
    
    done

done


for model_version in $model_version_list
    do

    python -u pretty_plots.py \
                --system mart \
                --system_name_in_csv_file mart_peptide_kd_averaged \
                --target_column'=-log10(Kd)' \
                --model hermes \
                --model_instance $model_version \
                --use_mt_structure $use_mt_structure \
                --show_wt_lines none
    
    python -u pretty_plots_respective_structures.py \
                --system mart \
                --system_name_in_csv_file 'mart_peptide_kd_{pdb}' \
                --target_column'=-log10(Kd)' \
                --model hermes \
                --model_instance $model_version \
                --use_mt_structure $use_mt_structure \
                --show_wt_lines none

done

