#!/bin/bash
#SBATCH --job-name=hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5__tcrdock
#SBATCH --account=stf
#SBATCH --partition=gpu-2080ti
#SBATCH --gres=gpu:1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=16:00:00
#SBATCH --mem=44G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=None
#SBATCH --export=all
#SBATCH -e /gscratch/spe/gvisan01/TCRdock-copy/slurm_logs/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5.err
#SBATCH -o /gscratch/spe/gvisan01/TCRdock-copy/slurm_logs/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5.out

python /gscratch/spe/gvisan01/TCRdock-copy/setup_for_alphafold.py                                 --targets_tsvfile /mmfs1/gscratch/spe/gvisan01/peptide_designs/ebv/hpvg_q5/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5.tsv                                 --output_dir /mmfs1/gscratch/spe/gvisan01/peptide_designs/ebv/hpvg_q5/hcnn_pyrosetta_annealing/tcrdock_output_noise=0p5                                 --new_docking

python -u /gscratch/spe/gvisan01/TCRdock-copy/run_prediction.py                                 --verbose                                 --targets /mmfs1/gscratch/spe/gvisan01/peptide_designs/ebv/hpvg_q5/hcnn_pyrosetta_annealing/tcrdock_output_noise=0p5/targets.tsv                                 --outfile_prefix /mmfs1/gscratch/spe/gvisan01/peptide_designs/ebv/hpvg_q5/hcnn_pyrosetta_annealing/tcrdock_output_noise=0p5/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5                                 --model_names model_2_ptm_ft4                                 --data_dir /gscratch/spe/gvisan01/TCRdock-copy/alphafold_params/                                 --model_params_files /gscratch/spe/gvisan01/TCRdock-copy/alphafold_params/params/tcrpmhc_run4_af_mhc_params_891.pkl

python /gscratch/spe/gvisan01/TCRdock-copy/add_pmhc_tcr_pae_to_tsvfile.py                                 --infile /mmfs1/gscratch/spe/gvisan01/peptide_designs/ebv/hpvg_q5/hcnn_pyrosetta_annealing/tcrdock_output_noise=0p5/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_final.tsv                                 --outfile /mmfs1/gscratch/spe/gvisan01/peptide_designs/ebv/hpvg_q5/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae.tsv                                --clobber