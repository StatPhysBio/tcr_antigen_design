#!/bin/bash
#SBATCH --job-name=magea3__sample_peptides_with_esmif
#SBATCH --account=spe
#SBATCH --partition=gpu-a40
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=23:00:00
#SBATCH --mem=44G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=gvisan01@uw.edu
#SBATCH --export=all
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/magea3__sample_peptides_with_esmif.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/magea3__sample_peptides_with_esmif.out

python sample_peptides_with_esmif.py
