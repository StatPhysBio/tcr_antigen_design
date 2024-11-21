

base_dir=/gscratch/spe/gvisan01/mutation_effect_predictions/hsiue_et_al/

model_instance_list='so3_convnet_base_finetuned_peptide_classification so3_convnet_noise=0p5_finetuned_peptide_classification' # 'so3_convnet_noise=0p5 so3_convnet_base_ensemble cgnet_base_ensemble cgnet_noise=0p5'
num_seq_per_target_list='10'

for model_instance in $model_instance_list
    do
    python -u ../plotting/saturation_mutagenesis_heatmap.py \
            --input_csv_file $base_dir'hcnn/zero_shot_predictions/'$model_instance'/hsiue_et_al_H2_sat_mut-'$model_instance'-use_mt_structure=0.csv' \
            --mutant_score_column log_proba_mt__minus__log_proba_wt \
            --show_mutants_in columns \
            --aa_order AILMFVPGWRHKDENQSTYC
done

# for model_instance in $model_instance_list
#     do
#     python -u ../plotting/saturation_mutagenesis_heatmap.py \
#             --input_csv_file $base_dir'hcnn/zero_shot_predictions/'$model_instance'/hsiue_et_al_H2_sat_mut-'$model_instance'-ensembles_of_pyrosetta_relaxed_structures.csv' \
#             --mutant_score_column log_proba_mt__minus__log_proba_wt \
#             --show_mutants_in columns \
#             --aa_order AILMFVPGWRHKDENQSTYC
# done

# for num_seq_per_target in $num_seq_per_target_list
#     do
#     python -u ../plotting/saturation_mutagenesis_heatmap.py \
#             --input_csv_file $base_dir'proteinmpnn/zero_shot_predictions/hsiue_et_al_H2_sat_mut-num_seq_per_target='$num_seq_per_target'-use_mt_structure=0.csv' \
#             --mutant_score_column log_p_mt__minus__log_p_wt \
#             --show_mutants_in columns \
#             --aa_order AILMFVPGWRHKDENQSTYC
# done


