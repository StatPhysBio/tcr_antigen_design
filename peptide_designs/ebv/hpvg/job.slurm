#!/bin/bash
#SBATCH --job-name=esmif_samples_0.7__tcrdock
#SBATCH --account=spe
#SBATCH --partition=gpu-a40
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=44G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=None
#SBATCH --export=all
#SBATCH -e /gscratch/spe/gvisan01/TCRdock-copy/slurm_logs/esmif_samples_0.7.err
#SBATCH -o /gscratch/spe/gvisan01/TCRdock-copy/slurm_logs/esmif_samples_0.7.out

python /gscratch/spe/gvisan01/TCRdock-copy/setup_for_alphafold.py                                 --targets_tsvfile /mmfs1/gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/esmif_0.7/esmif_samples_0.7.tsv                                 --output_dir /mmfs1/gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/esmif_0.7/tcrdock_output                                 --new_docking

python -u /gscratch/spe/gvisan01/TCRdock-copy/run_prediction.py                                 --verbose                                 --targets /mmfs1/gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/esmif_0.7/tcrdock_output/targets.tsv                                 --outfile_prefix /mmfs1/gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/esmif_0.7/tcrdock_output/esmif_samples_0.7                                 --model_names model_2_ptm_ft4                                 --data_dir /gscratch/spe/gvisan01/TCRdock-copy/alphafold_params/                                 --model_params_files /gscratch/spe/gvisan01/TCRdock-copy/alphafold_params/params/tcrpmhc_run4_af_mhc_params_891.pkl

python /gscratch/spe/gvisan01/TCRdock-copy/add_pmhc_tcr_pae_to_tsvfile.py                                 --infile /mmfs1/gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/esmif_0.7/tcrdock_output/esmif_samples_0.7_final.tsv                                 --outfile /mmfs1/gscratch/spe/gvisan01/tcr_pmhc/peptide_designs/ebv/hpvg/esmif_0.7/esmif_samples_0.7_w_pae.tsv                                --clobber