#!/bin/bash
#SBATCH --job-name=hsiue_hermes_py_050_6w51-filtered_HMTEVVRHC_4
#SBATCH --account=stf
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=6G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/hsiue_hermes_py_050_6w51-filtered_HMTEVVRHC_4.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/hsiue_hermes_py_050_6w51-filtered_HMTEVVRHC_4.out

python ../../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                                 --model_version hermes_py_050                                                 --pdb 6w51-filtered                                                 --pdbdir ./pdbs                                                 --chain F                                                 --peptide_resnum_start 1                                                 --sequence HMTEVVRHC                                                 --num_repeats 20                                                 --output_dir ./results/hermes_py_050/with_relaxation/                                                 --verbose 0                                                 --job 4