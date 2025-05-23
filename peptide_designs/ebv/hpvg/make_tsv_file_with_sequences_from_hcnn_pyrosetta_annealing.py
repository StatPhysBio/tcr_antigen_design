
import os
from glob import glob
import numpy as np
import pandas as pd


def get_true_sequences(sequences: np.ndarray, acceptance: np.ndarray):
    """Parse the accepted history from an annealing run"""
    true_seqs = np.zeros(
        shape=acceptance.shape, 
        dtype=sequences.dtype)
    
    true_seqs[0] = sequences[0]
    last_seq = sequences[0]
    for i,(seq,acc) in enumerate(zip(sequences[1:], acceptance[1:])):
        if acc:
            true_seqs[i+1] = seq
            last_seq = seq
        else:
            true_seqs[i+1] = last_seq
    return true_seqs


if __name__ == '__main__':

    for hcnn_model in ['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5']:

        run_files = glob(f'hcnn_pyrosetta_annealing/annealing_runs/*{hcnn_model}*.npy')

        print(f'Processing {len(run_files)} runs.')

        sequences = []
        pnlogps = []

        for run_file in run_files:

            run = np.load(run_file, allow_pickle=True)[()]

            try:
                true_sequences = get_true_sequences(run["sequences"][:-1], run["acceptance"])
                true_pnlogp = get_true_sequences(run['peptide pnlogp'], run['acceptance'])
            except KeyError:
                print(f"Skipping {run_file}")
                continue

            last_seq = true_sequences[-1].decode("utf-8")
            last_pnlogp = true_pnlogp[-1]
            
            sequences.append(last_seq)
            pnlogps.append(last_pnlogp)
        
        unique_sequences = []
        unique_pnlogps = []
        for seq, pnlogp in zip(sequences, pnlogps):
            if seq not in unique_sequences:
                unique_sequences.append(seq)
                unique_pnlogps.append(pnlogp)
        
        print(f'Unique sequences: {len(unique_sequences)}')

        ## now put these sequences in a tsv file for tcrdock to use
        header = 'organism	mhc_class	mhc	peptide	va	ja	cdr3a	vb	jb	cdr3b	pdbid	hcnn_model	pnlogp'
        template_row = 'human	1	B*35:01	{seq}	TRAV20*01	TRAJ58*01	CAVQDLGTSGSRLTF	TRBV9*01	TRBJ2-2*01	CASSARSGELFF	3mv7.pdb.human.MH1.B-35.A.C.DE	{hcnn_model}	{pnlogp}'
        
        with open(f'hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_{hcnn_model}.tsv', 'w+') as f:
            f.write(header + '\n')
            for seq, pnlogp in zip(sequences, pnlogps):
                f.write(template_row.format(seq=seq, hcnn_model=hcnn_model, pnlogp=pnlogp) + '\n')

        ## test the tsv file is read successfully
        df = pd.read_csv(f'hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_{hcnn_model}.tsv', sep='\t')
        print(df.head())

