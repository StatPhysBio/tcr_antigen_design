#!/bin/bash
#SBATCH --job-name=sample_peptides_from_blosum62__temperature=2.0__tcrdock
#SBATCH --account=spe
#SBATCH --partition=gpu-a40
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=16:00:00
#SBATCH --mem=96G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=None
#SBATCH --export=all
#SBATCH -e /gscratch/spe/gvisan01/TCRdock-copy/slurm_logs/sample_peptides_from_blosum62__temperature=2.0.err
#SBATCH -o /gscratch/spe/gvisan01/TCRdock-copy/slurm_logs/sample_peptides_from_blosum62__temperature=2.0.out

python /gscratch/spe/gvisan01/TCRdock-copy/setup_for_alphafold.py                                 --targets_tsvfile /mmfs1/gscratch/spe/gvisan01/peptide_designs/magea3_and_titin/titin/blosum62/sample_peptides_from_blosum62__temperature=2.0.tsv                                 --output_dir /mmfs1/gscratch/spe/gvisan01/peptide_designs/magea3_and_titin/titin/blosum62/tcrdock_output_T=2                                 --new_docking

python -u /gscratch/spe/gvisan01/TCRdock-copy/run_prediction.py                                 --verbose                                 --targets /mmfs1/gscratch/spe/gvisan01/peptide_designs/magea3_and_titin/titin/blosum62/tcrdock_output_T=2/targets.tsv                                 --outfile_prefix /mmfs1/gscratch/spe/gvisan01/peptide_designs/magea3_and_titin/titin/blosum62/tcrdock_output_T=2/sample_peptides_from_blosum62__temperature=2.0                                 --model_names model_2_ptm_ft4                                 --data_dir /gscratch/spe/gvisan01/TCRdock-copy/alphafold_params/                                 --model_params_files /gscratch/spe/gvisan01/TCRdock-copy/alphafold_params/params/tcrpmhc_run4_af_mhc_params_891.pkl

python /gscratch/spe/gvisan01/TCRdock-copy/add_pmhc_tcr_pae_to_tsvfile.py                                 --infile /mmfs1/gscratch/spe/gvisan01/peptide_designs/magea3_and_titin/titin/blosum62/tcrdock_output_T=2/sample_peptides_from_blosum62__temperature=2.0_final.tsv                                 --outfile /mmfs1/gscratch/spe/gvisan01/peptide_designs/magea3_and_titin/titin/blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae.tsv                                --clobber