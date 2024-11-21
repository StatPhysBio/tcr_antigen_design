
import os
import gzip, pickle
import numpy as np
import pandas as pd

from protein_holography_pytorch.utils.protein_naming import ol_to_ind_size, ind_to_ol_size

NUM_SAMPLES = 200


if __name__ == '__main__':

    os.makedirs('mhc_pwm', exist_ok=True)

    with gzip.open('../../all_structures/mhc_motif_pwms.pkl.gz', 'rb') as f:
        pwms = pickle.load(f)

    pwm = pwms['class_I']['B3501'][11]
    pwm = pwm / np.sum(pwm, axis=1, keepdims=True)

    samples = []
    for i in range(11):
        sampled_indices = np.random.multinomial(1, pwm[i, :], size=NUM_SAMPLES)
        characters = [ind_to_ol_size[np.argmax(sample)] for sample in sampled_indices]
        samples.append(characters)
    sequences = list(zip(*samples))
    sequences = [''.join(seq) for seq in sequences]
    sequences = list(set(sequences))

    ## now put these sequences in a tsv file for tcrdock to use
    header = 'organism	mhc_class	mhc	peptide	va	ja	cdr3a	vb	jb	cdr3b	pdbid	hcnn_model'
    template_row = 'human	1	B*35:01	{seq}	TRAV20*01	TRAJ58*01	CAVQDLGTSGSRLTF	TRBV9*01	TRBJ2-2*01	CASSARSGELFF	3mv7.pdb.human.MH1.B-35.A.C.DE	mhc_motif'
    
    with open('mhc_pwm/mhc_motif_peptides.tsv', 'w+') as f:
        f.write(header + '\n')
        for seq in sequences:
            f.write(template_row.format(seq=seq) + '\n')

    ## test the tsv file is read successfully
    df = pd.read_csv('mhc_pwm/mhc_motif_peptides.tsv', sep='\t')
    print(df.head())