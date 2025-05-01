#!/bin/bash
#SBATCH --job-name=mskcc_hermes_py_050_tcr1_NIVPMVATV_1
#SBATCH --account=stf
#SBATCH --partition=cpu-g2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=6G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_hermes_py_050_tcr1_NIVPMVATV_1.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_hermes_py_050_tcr1_NIVPMVATV_1.out

python ../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                                 --model_version hermes_py_050                                                 --pdb 5d2n-filtered                                                 --pdbdir ./pdbs                                                 --chain I                                                 --peptide_resnum_start 1                                                 --sequence NIVPMVATV                                                 --num_repeats 20                                                 --output_dir ./results/hermes_py_050/with_relaxation__runs_with_rosetta_scores/                                                 --verbose 0                                                 --job 1