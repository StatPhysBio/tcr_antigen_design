

proteinmpnn_dir=/gscratch/spe/gvisan01/ProteinMPNN-copy/mutation_effect_prediction/

pdbs='1ao7 1qse 1qsf'

base_dir='/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/tax/'

use_mt_structure='0'

model_version_list='v_48_030 v_48_020 v_48_002'

for pdb in $pdbs
    do

    echo $pdb

    for model_version in $model_version_list
        do

        python -u $proteinmpnn_dir'zero_shot_mutation_effect_prediction__faster.py' \
                        --csv_file $base_dir'tax_peptide_kd_'$pdb'.csv' \
                        --folder_with_pdbs $base_dir'pdbs/' \
                        --output_dir $base_dir'results/proteinmpnn_'$model_version \
                        --use_mt_structure $use_mt_structure \
                        --model_name $model_version \
                        --num_seq_per_target 10 \
                        --batch_size 10 \
                        --wt_pdb_column pdb \
                        --mutant_column mutants \
                        --mutant_chain_column chain \
                        --mutant_split_symbol '|'

        python -u pretty_plots.py \
                        --system tax \
                        --system_name_in_csv_file 'tax_peptide_kd_'$pdb \
                        --target_column'=neg_delta_log10_Kd' \
                        --model proteinmpnn \
                        --model_instance 'proteinmpnn_'$model_version \
                        --num_seq_per_target 10 \
                        --use_mt_structure $use_mt_structure \
                        --show_wt_lines both

    done
done


