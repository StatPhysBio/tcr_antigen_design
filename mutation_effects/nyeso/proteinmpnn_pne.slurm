#!/bin/bash
#SBATCH --job-name=mskcc_v_48_020
#SBATCH --account=spe
#SBATCH --partition=gpu-a40
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=23:00:00
#SBATCH --mem=48G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_v_48_020.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_v_48_020.out

bash proteinmpnn_pne.slurm

bash proteinmpnn_pne_tcr.slurm

bash proteinmpnn_pne_tcr_only.slurm
