#!/bin/bash
#SBATCH --job-name=mskcc_proteinmpnn
#SBATCH --account=spe
#SBATCH --partition=gpu-a40
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=23:00:00
#SBATCH --mem=48G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_proteinmpnn.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/mskcc_proteinmpnn.out

bash proteinmpnn_pne.sh

bash proteinmpnn_pne_tcr.sh

bash proteinmpnn_pne_tcr_only.sh
