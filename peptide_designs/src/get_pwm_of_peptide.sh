
model_version_list='hermes_py_000 hermes_py_050'

for model_version in $model_version_list
    do

    # ## nyeso

    # mkdir -p /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/nyeso/nyeso_full_copy/peptide_pwms

    # python get_pwm_of_peptide.py \
    #     --hermes_path /gscratch/spe/gvisan01/hermes/ \
    #     --model_version $model_version \
    #     --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/pdbs/2bnq.pdb \
    #     --chain C \
    #     --peptide_length 9 \
    #     --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/nyeso/nyeso_full_copy/peptide_pwms/pwm_2bnq_$model_version \
    #     --peptide_resnum_start 1

    # python get_pwm_of_peptide.py \
    #     --hermes_path /gscratch/spe/gvisan01/hermes/ \
    #     --model_version $model_version \
    #     --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/pdbs/2bnr.pdb \
    #     --chain C \
    #     --peptide_length 9 \
    #     --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/nyeso/nyeso_full_copy/peptide_pwms/pwm_2bnr_$model_version \
    #     --peptide_resnum_start 1

    # python get_pwm_of_peptide.py \
    #     --hermes_path /gscratch/spe/gvisan01/hermes/ \
    #     --model_version $model_version \
    #     --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/pdbs/fold_nyeso_no_template_model_1.pdb \
    #     --chain B \
    #     --peptide_length 9 \
    #     --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/nyeso/nyeso_full_copy/peptide_pwms/pwm_fold_nyeso_no_template_model_1_$model_version \
    #     --peptide_resnum_start 1

    # python get_pwm_of_peptide.py \
    #     --hermes_path /gscratch/spe/gvisan01/hermes/ \
    #     --model_version $model_version \
    #     --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/pdbs/fold_nyeso_yes_template_model_0.pdb \
    #     --chain B \
    #     --peptide_length 9 \
    #     --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/nyeso/nyeso_full_copy/peptide_pwms/pwm_fold_nyeso_yes_template_model_0_$model_version \
    #     --peptide_resnum_start 1

    # python get_pwm_of_peptide.py \
    #     --hermes_path /gscratch/spe/gvisan01/hermes/ \
    #     --model_version $model_version \
    #     --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/pdbs/fold_nyeso_9v_yes_template_model_0.pdb \
    #     --chain B \
    #     --peptide_length 9 \
    #     --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/nyeso/nyeso_full_copy/peptide_pwms/pwm_fold_nyeso_9v_yes_template_model_0_$model_version \
    #     --peptide_resnum_start 1

    # python get_pwm_of_peptide.py \
    #     --hermes_path /gscratch/spe/gvisan01/hermes/ \
    #     --model_version $model_version \
    #     --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/pdbs/fold_nyeso_9v_no_template_model_0.pdb \
    #     --chain B \
    #     --peptide_length 9 \
    #     --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/nyeso/nyeso_full_copy/peptide_pwms/pwm_fold_nyeso_9v_no_template_model_0_$model_version \
    #     --peptide_resnum_start 1
    

    ## ebv

    mkdir -p /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/peptide_pwms

    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/pdbs/3mv7.pdb.human.MH1.B-35.A.C.DE.pdb \
        --chain B \
        --peptide_length 11 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/peptide_pwms/pwm_3mv7_$model_version \
        --peptide_resnum_start 181
    
    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/pdbs/fold_ebv_no_template_model_0.pdb \
        --chain B \
        --peptide_length 11 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/peptide_pwms/pwm_fold_ebv_no_template_model_0_$model_version \
        --peptide_resnum_start 1

    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/pdbs/fold_ebv_yes_template_model_0.pdb \
        --chain B \
        --peptide_length 11 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/peptide_pwms/pwm_fold_ebv_yes_template_model_0_$model_version \
        --peptide_resnum_start 1

    mkdir -p /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg_q5/peptide_pwms

    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg_q5/pdbs/4prp.pdb.human.MH1.B-35.A.C.DE.pdb \
        --chain B \
        --peptide_length 11 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg_q5/peptide_pwms/pwm_4prp_$model_version \
        --peptide_resnum_start 181
    
    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg_q5/pdbs/fold_ebv_q5_no_template_model_0.pdb \
        --chain B \
        --peptide_length 11 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg_q5/peptide_pwms/pwm_fold_ebv_q5_no_template_model_0_$model_version \
        --peptide_resnum_start 1

    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg_q5/pdbs/fold_ebv_q5_yes_template_model_1.pdb \
        --chain B \
        --peptide_length 11 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg_q5/peptide_pwms/pwm_fold_ebv_q5_yes_template_model_1_$model_version \
        --peptide_resnum_start 1
    
    ## mage

    mkdir -p /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/magea3/peptide_pwms

    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/magea3/pdbs/5brz.pdb.human.MH1.A-01.A.C.DE.pdb \
        --chain B \
        --peptide_length 9 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/magea3/peptide_pwms/pwm_5brz_$model_version \
        --peptide_resnum_start 181
    
    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/magea3/pdbs/fold_mage_no_template_model_0.pdb \
        --chain B \
        --peptide_length 9 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/magea3/peptide_pwms/pwm_fold_mage_no_template_model_0_$model_version \
        --peptide_resnum_start 1
    
    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/magea3/pdbs/fold_mage_yes_template_seed_2_model_0.pdb \
        --chain B \
        --peptide_length 9 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/magea3/peptide_pwms/pwm_fold_mage_yes_template_seed_2_model_0_$model_version \
        --peptide_resnum_start 1

    mkdir -p /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/titin/peptide_pwms

    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/titin/pdbs/5bs0.pdb.human.MH1.A-01.A.C.DE.pdb \
        --chain B \
        --peptide_length 9 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/titin/peptide_pwms/pwm_5bs0_$model_version \
        --peptide_resnum_start 181
    
    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/titin/pdbs/fold_titin_no_template_model_0.pdb \
        --chain B \
        --peptide_length 9 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/titin/peptide_pwms/pwm_fold_titin_no_template_model_0_$model_version \
        --peptide_resnum_start 1
    
    python get_pwm_of_peptide.py \
        --hermes_path /gscratch/spe/gvisan01/hermes/ \
        --model_version $model_version \
        --pdbpath /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/titin/pdbs/fold_titin_yes_template_model_1.pdb \
        --chain B \
        --peptide_length 9 \
        --output_path_no_extension /gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/magea3_and_titin/titin/peptide_pwms/pwm_fold_titin_yes_template_model_1_$model_version \
        --peptide_resnum_start 1

done
