
import os

tcrdock_dir = '/gscratch/spe/gvisan01/TCRdock-copy'

## change the slurm setup according to the resouces available
os.makedirs('logs', exist_ok=True)
slurm_setup = "#!/bin/bash \n\
#SBATCH --job-name={system_identifier}__tcrdock \n\
#SBATCH --account=spe \n\
#SBATCH --partition=gpu-a40 \n\
#SBATCH --gres=gpu:1 \n\
#SBATCH --nodes=1 \n\
#SBATCH --ntasks-per-node=1 \n\
#SBATCH --time=12:00:00 \n\
#SBATCH --mem=64G \n\
#SBATCH --export=all \n\
#SBATCH -e logs/{system_identifier}__tcrdock.err \n\
#SBATCH -o logs/{system_identifier}__tcrdock.out"


output_dir = './'
inference_dir_template = os.path.join(output_dir, 'tcrdock_output_{temperature}')

template = 'sample_peptides_from_blosum62__temperature={temperature}'

for temperature in [1.0, 2.0, 3.0]:

    system_identifier = f'{template.format(temperature=temperature)}'
    tcrdock_input_file = f'{system_identifier}.tsv'
    inference_dir = inference_dir_template.format(temperature=temperature)
    os.makedirs(inference_dir, exist_ok=True)

    command_setup = f"python {tcrdock_dir}/setup_for_alphafold.py \
                                --targets_tsvfile {tcrdock_input_file} \
                                --output_dir {inference_dir} \
                                --new_docking"

    command_run = f"python -u {tcrdock_dir}/run_prediction.py \
                                --verbose \
                                --targets {os.path.join(inference_dir, 'targets.tsv')} \
                                --outfile_prefix {os.path.join(inference_dir, system_identifier)} \
                                --model_names model_2_ptm_ft4 \
                                --data_dir {tcrdock_dir}/alphafold_params/ \
                                --model_params_files {tcrdock_dir}/alphafold_params/params/tcrpmhc_run4_af_mhc_params_891.pkl"

    command_pae = f"python {tcrdock_dir}/add_pmhc_tcr_pae_to_tsvfile.py \
                                --infile {os.path.join(inference_dir, system_identifier + '_final.tsv')} \
                                --outfile {os.path.join(output_dir, system_identifier + '_w_pae.tsv')}\
                                --clobber"

    slurm_text = slurm_setup.format(system_identifier=system_identifier) + '\n\n' + command_setup + '\n\n' + command_run + '\n\n' + command_pae

    slurm_file = 'job.slurm'
    with open(slurm_file, 'w') as f:
        f.write(slurm_text)

    os.system(f"sbatch {slurm_file}")

