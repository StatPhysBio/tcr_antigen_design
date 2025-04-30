#!/bin/bash
#SBATCH --job-name=nyeso_hermes_py_050_2bnq_SLLMWIHQV_4
#SBATCH --account=stf
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=6G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/nyeso_hermes_py_050_2bnq_SLLMWIHQV_4.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/nyeso_hermes_py_050_2bnq_SLLMWIHQV_4.out

python ../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                                 --model_version hermes_py_050                                                 --pdb 2bnq                                                 --pdbdir ./pdbs                                                 --chain C                                                 --sequence SLLMWIHQV                                                 --num_repeats 20                                                 --output_dir ./results/hermes_py_050/with_relaxation__runs_with_rosetta_scores                                                 --verbose 0                                                 --job 4