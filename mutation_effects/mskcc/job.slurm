#!/bin/bash
#SBATCH --job-name=mskcc_hermes_py_050_tcr7_GRLKALCQR_4
#SBATCH --account=stf
#SBATCH --partition=gpu-l40s
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:30:00
#SBATCH --mem=5G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_hermes_py_050_tcr7_GRLKALCQR_4.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_hermes_py_050_tcr7_GRLKALCQR_4.out

python ../../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                                 --model_version hermes_py_050                                                 --pdb af3_tcr7                                                 --pdbdir ./pdbs                                                 --chain B                                                 --peptide_resnum_start 1                                                 --sequence GRLKALCQR                                                 --num_repeats 20                                                 --output_dir ./results/hermes_py_050/with_relaxation/                                                 --verbose 0                                                 --job 4