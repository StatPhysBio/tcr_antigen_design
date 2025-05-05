#!/bin/bash
#SBATCH --job-name=nyeso_hermes_py_050_2bnr_SLLMWITQV_4
#SBATCH --account=stf
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00
#SBATCH --mem=4G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/nyeso_hermes_py_050_2bnr_SLLMWITQV_4.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/nyeso_hermes_py_050_2bnr_SLLMWITQV_4.out

python ../../mutation_effects/src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                                 --model_version hermes_py_050                                                 --pdb 2bnr                                                 --pdbdir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/pdbs                                                 --chain C                                                 --sequence SLLMWITQV                                                 --num_repeats 20                                                 --output_dir ./hermes_scores/hermes_py_050/with_relaxation                                                 --verbose 0                                                 --job 4