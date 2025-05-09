#!/bin/bash
#SBATCH --job-name=nyeso_hermes_py_000_2bnq_GLFQWLAGV_4
#SBATCH --account=stf
#SBATCH --partition=gpu-l40s
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00
#SBATCH --mem=4G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/nyeso_hermes_py_000_2bnq_GLFQWLAGV_4.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/nyeso_hermes_py_000_2bnq_GLFQWLAGV_4.out

python ../../../mutation_effects/src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                                 --model_version hermes_py_000                                                 --pdb 2bnq                                                 --pdbdir /gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/pdbs                                                 --chain C                                                 --sequence GLFQWLAGV                                                 --peptide_resnum_start 1                                                 --num_repeats 20                                                 --output_dir ./hcnn_pyrosetta_annealing/hermes_py_000/with_relaxation                                                 --verbose 0                                                 --job 4