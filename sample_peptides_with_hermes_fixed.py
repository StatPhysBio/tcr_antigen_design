


import os
import numpy as np
import argparse
from typing import List

from hermes.inference import run_hermes_on_pdbfile_or_pyrosetta_pose

IDX_TO_AA = 'ACDEFGHIKLMNPQRSTVWY'
AA_TO_IDX = {aa: idx for idx, aa in enumerate(IDX_TO_AA)}


def sample_peptides_with_hermes_fixed(model_version: str, pdb_file: str, peptide_chain: str, temperature: float = 1.0, num_samples: int = 200, batch_size: int = 32) -> List[str]:

    ## run hermes
    df, _ = run_hermes_on_pdbfile_or_pyrosetta_pose(model_version, pdb_file, chain_and_sites_list=[peptide_chain], request=['logits'], batch_size=batch_size, ensemble_at_logits_level=True)

    len_peptide = len(df)
    num_aminoacids = len(IDX_TO_AA)

    ## extract logits
    logits = np.zeros((len_peptide, num_aminoacids))
    for i, row in df.iterrows():
        logits[i, :] = np.array([row[f'logit_{aa}'] for aa in IDX_TO_AA])
    
    ## apply temperature
    logits = logits / temperature

    ## convert to probabilities
    exp_logits = np.exp(logits - np.max(logits, axis=1, keepdims=True))
    pwm = exp_logits / np.sum(exp_logits, axis=1, keepdims=True)

    samples = []
    for i in range(len_peptide):
        sampled_indices = np.random.multinomial(1, pwm[i, :], size=num_samples)
        characters = [IDX_TO_AA[np.argmax(sample)] for sample in sampled_indices]
        samples.append(characters)
    
    sequences = list(zip(*samples))
    sequences = [''.join(seq) for seq in sequences]

    return sequences



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Sample peptides with HERMES fixed structure')

    parser.add_argument('-m', '--model_version', type=str, required=True,
                        help='Name of HERMES model you want to use. \
                              Use `hermes_py_000` or `hermes_bp_000` for the model trained without added noise to the atoms, \
                              `hermes_py_050` or `hermes_bp_050` for the model trained with 0.50 Angstrom noise to the atoms. \
                              Models with `_py_` in the name use pyrosetta to parse protein structures and were the ones used in the paper, \
                              whereas models with `_bp_` in the name use biopython.')
    
    parser.add_argument('-p', '--pdb_file', type=str, required=True,
                        help='Path to the PDB file containing the TCR-pMHC structure.')
    
    parser.add_argument('-c', '--peptide_chain', type=str, required=True,
                        help='Chain identifier for the peptide in the PDB file.')

    parser.add_argument('-o', '--output_file', type=str, required=True,
                        help='Path to the output .txt file where sampled peptide sequences will be saved (one sequence per line).')

    parser.add_argument('-t', '--temperature', type=float, default=1.0,
                            help='Temperature for sampling from the probability distribution (default: 1.0, lower is more conservative, higher is more diverse).')
    
    parser.add_argument('-n', '--num_samples', type=int, default=200,
                         help='Number of peptide samples to generate (note that they are *not* guaranteed to be unique). Default: 200.')
    
    parser.add_argument('-b', '--batch_size', type=int, default=32,
                         help='Batch size (number of sites) for running HERMES inference. Default: 32, no need to change this except for handling memory issues.')
    
    args = parser.parse_args()

    sequences = sample_peptides_with_hermes_fixed(args.model_version, args.pdb_file, args.peptide_chain, args.temperature, args.num_samples, args.batch_size)
    with open(args.output_file, 'w') as f:
        for seq in sequences:
            f.write(seq + '\n')
    print(f'Sampled {len(sequences)} peptide sequences and saved them to {args.output_file}.')

