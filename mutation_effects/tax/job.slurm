#!/bin/bash
#SBATCH --job-name=tax_hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi_1qsf_LLFGHPVAV_4
#SBATCH --account=spe
#SBATCH --partition=ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=01:00:00
#SBATCH --mem=4G
#SBATCH -e ./__slurm_logs/tax_hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi_1qsf_LLFGHPVAV_4.err
#SBATCH -o ./__slurm_logs/tax_hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi_1qsf_LLFGHPVAV_4.out

python ../../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                                 --model_version hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi                                                 --pdb 1qsf                                                 --pdbdir ./pdbs                                                 --chain C                                                 --sequence LLFGHPVAV                                                 --num_repeats 20                                                 --output_dir ./results/hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi/with_relaxation/                                                 --verbose 0                                                 --job 4