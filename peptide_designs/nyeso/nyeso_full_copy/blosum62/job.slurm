#!/bin/bash 
#SBATCH --job-name=sample_peptides_from_blosum62__temperature=3.0__tcrdock 
#SBATCH --account=spe 
#SBATCH --partition=gpu-a40 
#SBATCH --gres=gpu:1 
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=1 
#SBATCH --time=12:00:00 
#SBATCH --mem=64G 
#SBATCH --mail-type=ALL 
#SBATCH --mail-user=gvisan01@uw.edu 
#SBATCH --export=all 
#SBATCH -e logs/sample_peptides_from_blosum62__temperature=3.0__tcrdock.err 
#SBATCH -o logs/sample_peptides_from_blosum62__temperature=3.0__tcrdock.out

python /gscratch/spe/gvisan01/TCRdock-copy/setup_for_alphafold.py                                 --targets_tsvfile sample_peptides_from_blosum62__temperature=3.0.tsv                                 --output_dir ./tcrdock_output_3.0                                 --new_docking

python -u /gscratch/spe/gvisan01/TCRdock-copy/run_prediction.py                                 --verbose                                 --targets ./tcrdock_output_3.0/targets.tsv                                 --outfile_prefix ./tcrdock_output_3.0/sample_peptides_from_blosum62__temperature=3.0                                 --model_names model_2_ptm_ft4                                 --data_dir /gscratch/spe/gvisan01/TCRdock-copy/alphafold_params/                                 --model_params_files /gscratch/spe/gvisan01/TCRdock-copy/alphafold_params/params/tcrpmhc_run4_af_mhc_params_891.pkl

python /gscratch/spe/gvisan01/TCRdock-copy/add_pmhc_tcr_pae_to_tsvfile.py                                 --infile ./tcrdock_output_3.0/sample_peptides_from_blosum62__temperature=3.0_final.tsv                                 --outfile ./sample_peptides_from_blosum62__temperature=3.0_w_pae.tsv                                --clobber