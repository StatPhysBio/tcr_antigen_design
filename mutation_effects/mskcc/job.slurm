#!/bin/bash
#SBATCH --job-name=mskcc_hermes_py_050_tcr2_NLVPMVATV_4
#SBATCH --account=stf
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:00:00
#SBATCH --mem=6G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_hermes_py_050_tcr2_NLVPMVATV_4.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_hermes_py_050_tcr2_NLVPMVATV_4.out

python ../../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                                 --model_version hermes_py_050                                                 --pdb TCR2_T00000_A0201_NLVPMVATV_0_model_1_model_2_ptm_ft4                                                 --pdbdir ./pdbs                                                 --chain A                                                 --peptide_resnum_start 376                                                 --sequence NLVPMVATV                                                 --num_repeats 20                                                 --output_dir ./results/hermes_py_050/with_relaxation/                                                 --verbose 0                                                 --job 4