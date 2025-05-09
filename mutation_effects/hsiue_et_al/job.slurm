#!/bin/bash
#SBATCH --job-name=hsiue_et_al_H2_sat_mut_af3_yes_template.csv_hermes_py_050_fold_hsiue_yes_template_model_0_HMTEVVRHC_4
#SBATCH --account=stf
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/hsiue_et_al_H2_sat_mut_af3_yes_template.csv_hermes_py_050_fold_hsiue_yes_template_model_0_HMTEVVRHC_4.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/hsiue_et_al_H2_sat_mut_af3_yes_template.csv_hermes_py_050_fold_hsiue_yes_template_model_0_HMTEVVRHC_4.out

python ../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                             --model_version hermes_py_050                                             --pdb fold_hsiue_yes_template_model_0                                             --pdbdir ./pdbs                                             --chain C                                             --peptide_resnum_start 1                                             --sequence HMTEVVRHC                                             --num_repeats 20                                             --output_dir ./results/hermes_py_050/with_relaxation/                                             --verbose 0                                             --job 4