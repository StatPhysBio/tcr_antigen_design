
import os, sys
import numpy as np
import pandas as pd

import logomaker

# get numpy largest integer
MIN = -np.iinfo(np.int32).max

np.random.seed(42)

from constants import WT_SEQ


from Bio.Align import substitution_matrices

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

BLOSUM62 = substitution_matrices.load('BLOSUM62')

BLOSUM62_DICT = {aa1: {aa2: BLOSUM62[aa1, aa2] for aa2 in AMINO_ACIDS} for aa1 in AMINO_ACIDS}

BLOSUM62_NO_SELF_DICT = {aa1: {aa2: BLOSUM62_DICT[aa1][aa2] if aa1 != aa2 else MIN for aa2 in AMINO_ACIDS} for aa1 in AMINO_ACIDS}

# compute softmax of blosum62
SOFTMAX_BLOSUM62_DICT = {aa1: {aa2: np.exp(BLOSUM62_DICT[aa1][aa2]) for aa2 in AMINO_ACIDS} for aa1 in AMINO_ACIDS}
SOFTMAX_BLOSUM62_DICT = {aa1: {aa2: SOFTMAX_BLOSUM62_DICT[aa1][aa2] / sum(list(SOFTMAX_BLOSUM62_DICT[aa1].values())) for aa2 in AMINO_ACIDS} for aa1 in AMINO_ACIDS}

SOFTMAX_BLOSUM62_NO_SELF_DICT = {aa1: {aa2: np.exp(BLOSUM62_NO_SELF_DICT[aa1][aa2]) for aa2 in AMINO_ACIDS} for aa1 in AMINO_ACIDS}
SOFTMAX_BLOSUM62_NO_SELF_DICT = {aa1: {aa2: SOFTMAX_BLOSUM62_NO_SELF_DICT[aa1][aa2] / sum(list(SOFTMAX_BLOSUM62_NO_SELF_DICT[aa1].values())) for aa2 in AMINO_ACIDS} for aa1 in AMINO_ACIDS}


if __name__ == '__main__':

    os.makedirs('blosum62', exist_ok=True)

    NUM_TARGET_SAMPLES = 200

    header = "organism	mhc_class	mhc	peptide	va	ja	cdr3a	vb	jb	cdr3b	pdbid	hcnn_model"
    template_row = "human	1	A*01:01	{pep_seq}	TRAV21*01	TRAJ28*01	CAVRPGGAGPFFVVF	TRBV5-1*01	TRBJ2-7*01	CASSFNMATGQYF	5bs0.pdb.human.MH1.A-01.A.C.DE	blosum62"

    for temperature in [1.0, 2.0]:

        softmax_blosum62_dict = {aa1: {aa2: np.exp(BLOSUM62_DICT[aa1][aa2] / temperature) for aa2 in AMINO_ACIDS} for aa1 in AMINO_ACIDS}
        softmax_blosum62_dict = {aa1: {aa2: softmax_blosum62_dict[aa1][aa2] / sum(list(softmax_blosum62_dict[aa1].values())) for aa2 in AMINO_ACIDS} for aa1 in AMINO_ACIDS}

        pep_sequences_unique = set()

        while len(pep_sequences_unique) < NUM_TARGET_SAMPLES:
            # sample from blosum

            pep_seq = []

            for j in range(len(WT_SEQ)):
                pep_seq.append(np.random.choice(AMINO_ACIDS, p=list(softmax_blosum62_dict[WT_SEQ[j]].values())))

            pep_seq = ''.join(pep_seq)

            pep_sequences_unique.add(pep_seq)

        pep_sequences_unique = list(pep_sequences_unique)
        
        # make logo of sequences

        # make pwm
        pwm = np.zeros((len(pep_sequences_unique[0]), len(AMINO_ACIDS)))

        for i, pep_seq in enumerate(pep_sequences_unique):
            for j, aa in enumerate(pep_seq):
                pwm[j, AMINO_ACIDS.index(aa)] += 1
        
        pwm = pwm / len(pep_sequences_unique)
        pwm_df = pd.DataFrame(pwm, columns=AMINO_ACIDS)

        # make logo
        logo = logomaker.Logo(pwm_df)
        logo.ax.set_title('Sampled Peptides from BLOSUM62')
        logo.fig.savefig(f'blosum62/sample_peptides_from_blosum62__temperature={temperature}.png')


        # make tsv file
        with open(f'blosum62/sample_peptides_from_blosum62__temperature={temperature}.tsv', 'w+') as f:
            f.write(header + '\n')
            for pep_seq in pep_sequences_unique:
                f.write(template_row.format(pep_seq=pep_seq) + '\n')

        # make sure it reads correctly
        df = pd.read_csv(f'blosum62/sample_peptides_from_blosum62__temperature={temperature}.tsv', sep='\t')
        print(df.head())






