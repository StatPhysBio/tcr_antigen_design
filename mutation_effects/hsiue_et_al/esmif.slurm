#!/bin/bash
#SBATCH --job-name=hsiue_esmif
#SBATCH --account=spe
#SBATCH --partition=gpu-a40
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=64G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/hsiue_esmif.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/hsiue_esmif.out

base_dir='/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/hsiue_et_al/'

python -u ../src/score_peptides_with_esmif.py \
            --pdbdir $base_dir'pdbs' \
            --output_csv_filepath $base_dir'results/esmif/hsiue_et_al_H2_sat_mut-esmif-use_mt_structure=0.csv' \
            --csv_filepath $base_dir'hsiue_et_al_H2_sat_mut.csv' \
            --pdb_column wt_pdb \
            --chain_column mutant_chain \
            --peptide_column sequence
