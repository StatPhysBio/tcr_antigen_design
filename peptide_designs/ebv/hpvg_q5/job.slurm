#!/bin/bash
#SBATCH --job-name=ebv_hermes_py_050_4prp.pdb.human.MH1.B-35.A.C.DE_YPFGPFDYLLY_4
#SBATCH --account=stf
#SBATCH --partition=cpu-g2-mem2x
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=03:00:00
#SBATCH --mem=5G
#SBATCH -e /gscratch/scrubbed/gvisan01/slurm_logs/ebv_hermes_py_050_4prp.pdb.human.MH1.B-35.A.C.DE_YPFGPFDYLLY_4.err
#SBATCH -o /gscratch/scrubbed/gvisan01/slurm_logs/ebv_hermes_py_050_4prp.pdb.human.MH1.B-35.A.C.DE_YPFGPFDYLLY_4.out

python ../../../mutation_effects/src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py                                                 --model_version hermes_py_050                                                 --pdb 4prp.pdb.human.MH1.B-35.A.C.DE                                                 --pdbdir ./pdbs                                                 --chain B                                                 --sequence YPFGPFDYLLY                                                 --peptide_resnum_start 181                                                 --num_repeats 20                                                 --output_dir ./hcnn_pyrosetta_annealing/hermes_py_050/with_relaxation                                                 --verbose 0                                                 --job 4