


import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO

import logomaker


NUM_SAMPLES = 200


if __name__ == '__main__':

    curr_dir = os.path.dirname(os.path.realpath(__file__))

    folder_with_pdbs = os.path.join(curr_dir, 'pdbs')
    output_dir_template = os.path.join(curr_dir, 'proteinmpnn')

    model_names = ['v_48_020', 'v_48_002']

    for sampling_temp in [0.1, 0.4, 0.7]:

        for model_name in model_names:

            command = f"python -u /gscratch/spe/gvisan01/ProteinMPNN-copy/peptide_scoring_and_design/design_peptides.py \
                                        --folder_with_pdbs {folder_with_pdbs} \
                                        --output_dir_template {output_dir_template} \
                                        --sampling_temp {sampling_temp} \
                                        --num_seq_per_target {NUM_SAMPLES} \
                                        --seed 37 \
                                        --batch_size 25 \
                                        --peptide_chain B \
                                        --peptide_resnum_start 1 \
                                        --peptide_length 11 \
                                        --model_name {model_name}"

            os.system(command)

            output_dir = output_dir_template + f'_{model_name}' + f'_{sampling_temp}'

            print(model_name)
            pdbids = [f[:-4] for f in os.listdir(folder_with_pdbs) if f.endswith('.pdb')]
            sequences_per_pdbid = {}
            for pdbid in pdbids:
                print(pdbid)
                sequences = []
                with open(os.path.join(output_dir, 'seqs', f'{pdbid}.fa'), 'r') as f:
                    for i, record in enumerate(SeqIO.parse(f, 'fasta')):
                        if i > 0:
                            sequences.append(str(record.seq))
                
                sequences = list(set(sequences))
                print(f'Total unique sequences: {len(sequences)}')
                sequences_per_pdbid[pdbid] = sequences
                    
                # make pwm, save it as image
                pwm = logomaker.alignment_to_matrix(sequences)
                pwm_df = pd.DataFrame(pwm, columns=list('ACDEFGHIKLMNPQRSTVWY'))
                for i in range(pwm_df.shape[0]):
                    pwm_df.iloc[i] = pwm_df.iloc[i] / pwm_df.iloc[i].sum()

                # set nan values to zero
                pwm_df = pwm_df.fillna(0)
                
                fig, ax = plt.subplots(figsize=(10, 2))
                logomaker.Logo(pwm_df, ax=ax)
                plt.tight_layout()
                plt.savefig(os.path.join(output_dir, 'seqs', f'{pdbid}_logo.png'))
                plt.close()

            print()

            ## now put these sequences in a tsv file for tcrdock to use
            header = 'organism	mhc_class	mhc	peptide	va	ja	cdr3a	vb	jb	cdr3b	pdbid	hcnn_model'
            template_row = 'human	1	B*35:01	{seq}	TRAV20*01	TRAJ58*01	CAVQDLGTSGSRLTF	TRBV9*01	TRBJ2-2*01	CASSARSGELFF	{pdbid}	proteinmpnn_{model_name}'
            
            with open(os.path.join(output_dir, f'proteinmpnn_samples_{model_name}_{sampling_temp}.tsv'), 'w+') as f:
                f.write(header + '\n')
                written_sequences = set()
                for pdbid, sequences in sequences_per_pdbid.items():
                    for seq in sequences:
                        if seq in written_sequences:
                            continue
                        else:
                            f.write(template_row.format(seq=seq, pdbid=pdbid, model_name=model_name) + '\n')
                            written_sequences.add(seq)

            ## test the tsv file is read successfully
            df = pd.read_csv(os.path.join(output_dir, f'proteinmpnn_samples_{model_name}_{sampling_temp}.tsv'), sep='\t')
            print(df.head())
