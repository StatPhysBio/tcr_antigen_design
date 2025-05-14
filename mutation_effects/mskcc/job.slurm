#!/bin/bash
#SBATCH --job-name=mskcc_tcr2_ec50_sat_mut_af3_no_template_hermes_py_050_NLVPMVATR_4
#SBATCH --account=stf
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00
#SBATCH --mem=5G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_tcr2_ec50_sat_mut_af3_no_template_hermes_py_050_NLVPMVATR_4.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_tcr2_ec50_sat_mut_af3_no_template_hermes_py_050_NLVPMVATR_4.out

python ../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                             --model_version hermes_py_050                                             --pdb fold_tcr2_no_template_model_1                                             --pdbdir ./pdbs                                             --chain B                                             --peptide_resnum_start 1                                             --sequence NLVPMVATR                                             --num_repeats 20                                             --output_dir ./results/hermes_py_050/with_relaxation/                                             --verbose 0                                             --job 4