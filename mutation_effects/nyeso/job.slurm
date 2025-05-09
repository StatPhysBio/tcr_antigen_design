#!/bin/bash
#SBATCH --job-name=nyeso_hermes_py_050_fold_nyeso_9v_no_template_model_0_SLLMWIHQV_4
#SBATCH --account=stf
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=4G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/nyeso_hermes_py_050_fold_nyeso_9v_no_template_model_0_SLLMWIHQV_4.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/nyeso_hermes_py_050_fold_nyeso_9v_no_template_model_0_SLLMWIHQV_4.out

python ../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                             --model_version hermes_py_050                                             --pdb fold_nyeso_9v_no_template_model_0                                             --pdbdir ./pdbs                                             --chain B                                             --sequence SLLMWIHQV                                             --num_repeats 20                                             --output_dir ./results/hermes_py_050/with_relaxation                                             --verbose 0                                             --job 4