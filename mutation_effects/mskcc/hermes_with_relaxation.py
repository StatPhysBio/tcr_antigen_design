
import os
import numpy as np
import pandas as pd
import argparse

from constants import TCR_TO_PDB, PDB_TO_PEP_INFO

# done TCRS: 1, 2, 3, 

TCRS = ['4']


SLURM_SETUP = "#!/bin/bash\n\
#SBATCH --job-name={system_identifier}\n\
#SBATCH --account={account}\n\
#SBATCH --partition={partition}\n{gpu_text}\
#SBATCH --nodes=1\n\
#SBATCH --ntasks-per-node={num_cores}\n\
#SBATCH --time={walltime}\n\
#SBATCH --mem={memory}\n{email_text}\
#SBATCH -e {errfile}\n\
#SBATCH -o {outfile}"


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-m',  '--model_version', type=str, nargs='+', default=['hermes_py_000', 'hermes_py_050']) #, 'hermes_py_000_ft_skempi_no_tcrpmhc_ddg_bi', 'hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi']) # hermes_py_000 hermes_py_050 hermes_py_000_ft_skempi_no_tcrpmhc_ddg_bi hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi')
    parser.add_argument('--num_repeats', type=int, default=20)
    parser.add_argument('--total_relaxations', type=int, default=100)

    parser.add_argument('-A',  '--account', type=str, default='stf')
    parser.add_argument('-P',  '--partition', type=str, default='gpu-l40')
    parser.add_argument('-G',  '--use_gpu', type=int, default=0, choices=[0, 1])
    parser.add_argument('-C',  '--num_cores', type=int, default=1)
    parser.add_argument('-W',  '--walltime', type=str, default='02:00:00')
    parser.add_argument('-M',  '--memory', type=str, default='8G')
    parser.add_argument('-E',  '--send_emails', type=int, default=0, choices=[0, 1])
    parser.add_argument('-EA', '--email_address', type=str, default='gvisan01@uw.edu')

    args = parser.parse_args()

    logs_path = '/gscratch/scrubbed/gvisan01/slurm_logs'
    os.makedirs(logs_path, exist_ok=True)
    
    if args.use_gpu:
        gpu_text = '#SBATCH --gres=gpu:1\n'
    else:
        gpu_text = ''
    
    if args.send_emails:
        email_text = f'#SBATCH --mail-type=ALL\n#SBATCH --mail-user={args.email_address}\n#SBATCH --export=all\n'
    else:
        email_text = ''
    
    for model_version in args.model_version:
        os.makedirs(f'./results/{model_version}/with_relaxation', exist_ok=True)
    
        for tcr in TCRS:
            pdb = TCR_TO_PDB[tcr]
            df = pd.read_csv(f'mskcc_tcr{tcr}_ec50_sat_mut.csv')
            for i, row in df.iterrows():
                sequence = row['sequence']
                num_jobs = args.total_relaxations // args.num_repeats

                for j in range(num_jobs):

                    identifier = f'mskcc_{model_version}_tcr{tcr}_{sequence}_{j}'

                    slurm_text = SLURM_SETUP.format(system_identifier=identifier,
                                                    account=args.account,
                                                    partition=args.partition,
                                                    gpu_text=gpu_text,
                                                    num_cores=args.num_cores,
                                                    walltime=args.walltime,
                                                    memory=args.memory,
                                                    email_text=email_text,
                                                    errfile=os.path.join(logs_path, f"{identifier}.err"),
                                                    outfile=os.path.join(logs_path, f"{identifier}.out"))
                    slurm_text += '\n\n'

                    slurm_text += f"python ../../src/get_score_of_requested_peptide_in_pdb_of_tcr_pmhc_with_relaxation.py \
                                                --model_version {model_version} \
                                                --pdb {pdb} \
                                                --pdbdir ./pdbs \
                                                --chain {PDB_TO_PEP_INFO[pdb][0]} \
                                                --peptide_resnum_start {PDB_TO_PEP_INFO[pdb][1]} \
                                                --sequence {sequence} \
                                                --num_repeats {args.num_repeats} \
                                                --output_dir ./results/{model_version}/with_relaxation/ \
                                                --verbose 0 \
                                                --job {j}"


                    slurm_file = 'job.slurm'
                    with open(slurm_file, 'w') as f:
                        f.write(slurm_text)
                    
                    os.system(f'sbatch {slurm_file}')

