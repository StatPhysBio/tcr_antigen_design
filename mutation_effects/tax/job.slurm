#!/bin/bash
#SBATCH --job-name=tax_hermes_py_050_1qsf_LLFGHPVAV_4
#SBATCH --account=stf
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=6G
#SBATCH -e ./__slurm_logs/tax_hermes_py_050_1qsf_LLFGHPVAV_4.err
#SBATCH -o ./__slurm_logs/tax_hermes_py_050_1qsf_LLFGHPVAV_4.out

python ../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                                 --model_version hermes_py_050                                                 --pdb 1qsf                                                 --pdbdir ./pdbs                                                 --chain C                                                 --sequence LLFGHPVAV                                                 --num_repeats 20                                                 --output_dir ./results/hermes_py_050/with_relaxation__runs_with_rosetta_scores/                                                 --verbose 0                                                 --job 4