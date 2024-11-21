

proteinmpnn_dir=/gscratch/spe/gvisan01/ProteinMPNN-copy/mutation_effect_prediction/

base_dir='/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/chandran_et_al/'

tcrs='3 4'

use_mt_structure='0'

model_version_list='v_48_030 v_48_020 v_48_002'

for tcr in $tcrs
    do

    for model_version in $model_version_list
        do

        python -u $proteinmpnn_dir'zero_shot_mutation_effect_prediction__faster.py' \
                        --csv_file $base_dir'chandran_et_al_peptide_A_and_G_scans_TCR'$tcr'.csv' \
                        --folder_with_pdbs $base_dir'pdbs/' \
                        --output_dir $base_dir'results/proteinmpnn_'$model_version \
                        --use_mt_structure $use_mt_structure \
                        --model_name $model_version \
                        --num_seq_per_target 10 \
                        --batch_size 10 \
                        --wt_pdb_column wt_pdb \
                        --mutant_column mutant \
                        --mutant_chain_column mutant_chain

        python -u pretty_plots.py \
                    --system chandran_et_al \
                    --system_name_in_csv_file 'chandran_et_al_peptide_A_and_G_scans_TCR'$tcr \
                    --target_column'=TFN_alpha_mt/TFN_alpha_wt' \
                    --model proteinmpnn \
                    --model_instance 'proteinmpnn_'$model_version \
                    --num_seq_per_target 10 \
                    --use_mt_structure $use_mt_structure \
                    --show_wt_lines both \
                    --wt_value_target 100.0 \
                    --wt_value_pred 0.0

    done

done

