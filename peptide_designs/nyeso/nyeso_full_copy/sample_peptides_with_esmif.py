
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import SeqIO

import logomaker

ESM_PATH = '/gscratch/spe/gvisan01/esm/esm'

NUM_SAMPLES = 100


if __name__ == '__main__':

    curr_dir = os.path.dirname(os.path.realpath(__file__))

    pdbdir = '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/pdbs/'
    pdbids = ['2bnq', '2bnr']
    output_dir_template = os.path.join(curr_dir, 'esmif')


    for sampling_temp in [0.1, 0.4, 0.7]:
            
        output_dir = output_dir_template + f'_{sampling_temp}'
        os.makedirs(output_dir, exist_ok=True)

        output_fasta_path = os.path.join(output_dir, 'seqs', 'seqs.fasta')

        outfile_handle = open(os.path.join(output_dir, f'esmif_samples_{sampling_temp}.tsv'), 'w+')

        for pdb_i, pdbid in enumerate(pdbids):

            pdbpath = os.path.join(pdbdir, pdbid + '.pdb')

            command = f'python {ESM_PATH}/examples/inverse_folding/sample_sequences.py \
                                    {pdbpath} \
                                    --chain C \
                                    --temperature {sampling_temp} \
                                    --num-samples {NUM_SAMPLES} \
                                    --outpath {output_fasta_path} \
                                    --multichain-backbone'

            os.system(command)

            sequences = []
            with open(output_fasta_path, 'r') as f:
                for i, record in enumerate(SeqIO.parse(f, 'fasta')):
                    if i > 0:
                        sequences.append(str(record.seq))
            
            sequences = list(set(sequences))
            print(f'Total unique sequences: {len(sequences)}')
            
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
            plt.savefig(os.path.join(output_dir, 'seqs', f'seqs_logo_{pdbid}.png'))
            plt.close()

            print()

            ## now put these sequences in a tsv file for tcrdock to use
            header = 'organism	mhc_class	mhc	peptide	va	ja	cdr3a	vb	jb	cdr3b	pdbid	hcnn_model'
            template_row = 'human	1	A*02:01	{seq}	TRAV21*01	TRAJ6*01	CAVRPTSGGSYIPTF	TRBV6-5*01	TRBJ2-2*01	CASSYVGNTGELFF	{pdbid}	esmif'
            
            if pdb_i == 0:
                outfile_handle.write(header + '\n')
            
            written_sequences = set()
            for seq in sequences:
                if seq in written_sequences:
                    continue
                else:
                    outfile_handle.write(template_row.format(seq=seq, pdbid=pdbid) + '\n')
                    written_sequences.add(seq)
        
        outfile_handle.close()

        ## test the tsv file is read successfully
        df = pd.read_csv(os.path.join(output_dir, f'esmif_samples_{sampling_temp}.tsv'), sep='\t')
        print(df.head())


