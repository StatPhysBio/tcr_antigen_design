
import os, sys
import numpy as np
import pandas as pd

import logomaker

from constants import WT_SEQ

# get numpy largest integer
MIN = -np.iinfo(np.int32).max

np.random.seed(42)


from Bio.Align import substitution_matrices

AMINO_ACIDS = list('ACDEFGHIKLMNPQRSTVWY')

BLOSUM62_DICT = {
    'C': {'C': 9, 'S': -1, 'T': -1, 'A': 0, 'G': -3, 'P': -3, 'D': -3, 'E': -4, 'Q': -3, 'N': -3, 'H': -3, 'R': -3, 'K': -3, 'M': -1, 'I': -1, 'L': -1, 'V': -1, 'W': -2, 'Y': -2, 'F': -2},
    'S': {'C': -1, 'S': 4, 'T': 1, 'A': 1, 'G': 0, 'P': -1, 'D': 0, 'E': 0, 'Q': 0, 'N': 1, 'H': -1, 'R': -1, 'K': 0, 'M': -1, 'I': -2, 'L': -2, 'V': -2, 'W': -3, 'Y': -2, 'F': -2},
    'T': {'C': -1, 'S': 1, 'T': 5, 'A': 0, 'G': -2, 'P': -1, 'D': -1, 'E': -1, 'Q': -1, 'N': 0, 'H': -2, 'R': -1, 'K': -1, 'M': -1, 'I': -1, 'L': -1, 'V': 0, 'W': -2, 'Y': -2, 'F': -2},
    'A': {'C': 0, 'S': 1, 'T': 0, 'A': 4, 'G': 0, 'P': -1, 'D': -2, 'E': -1, 'Q': -1, 'N': -2, 'H': -2, 'R': -1, 'K': -1, 'M': -1, 'I': -1, 'L': -1, 'V': 0, 'W': -3, 'Y': -2, 'F': -2},
    'G': {'C': -3, 'S': 0, 'T': -2, 'A': 0, 'G': 6, 'P': -2, 'D': -1, 'E': -2, 'Q': -2, 'N': 0, 'H': -2, 'R': -2, 'K': -2, 'M': -3, 'I': -4, 'L': -4, 'V': -3, 'W': -2, 'Y': -3, 'F': -3},
    'P': {'C': -3, 'S': -1, 'T': -1, 'A': -1, 'G': -2, 'P': 7, 'D': -1, 'E': -1, 'Q': -1, 'N': -1, 'H': -2, 'R': -2, 'K': -1, 'M': -2, 'I': -3, 'L': -3, 'V': -2, 'W': -4, 'Y': -3, 'F': -4},
    'D': {'C': -3, 'S': 0, 'T': -1, 'A': -2, 'G': -1, 'P': -1, 'D': 6, 'E': 2, 'Q': 0, 'N': 1, 'H': -1, 'R': -2, 'K': -1, 'M': -3, 'I': -3, 'L': -4, 'V': -3, 'W': -4, 'Y': -3, 'F': -3},
    'E': {'C': -4, 'S': 0, 'T': -1, 'A': -1, 'G': -2, 'P': -1, 'D': 2, 'E': 5, 'Q': 2, 'N': 0, 'H': 0, 'R': 0, 'K': 1, 'M': -2, 'I': -3, 'L': -3, 'V': -2, 'W': -3, 'Y': -2, 'F': -3},
    'Q': {'C': -3, 'S': 0, 'T': -1, 'A': -1, 'G': -2, 'P': -1, 'D': 0, 'E': 2, 'Q': 5, 'N': 0, 'H': 0, 'R': 1, 'K': 1, 'M': 0, 'I': -3, 'L': -2, 'V': -2, 'W': -2, 'Y': -1, 'F': -3},
    'N': {'C': -3, 'S': 1, 'T': 0, 'A': -2, 'G': 0, 'P': -2, 'D': 1, 'E': 0, 'Q': 0, 'N': 6, 'H': 1, 'R': 0, 'K': 0, 'M': -2, 'I': -3, 'L': -3, 'V': -3, 'W': -4, 'Y': -2, 'F': -3},
    'H': {'C': -3, 'S': -1, 'T': -2, 'A': -2, 'G': -2, 'P': -2, 'D': -1, 'E': 0, 'Q': 0, 'N': 1, 'H': 8, 'R': 0, 'K': -1, 'M': -2, 'I': -3, 'L': -3, 'V': -3, 'W': -2, 'Y': 2, 'F': -1},
    'R': {'C': -3, 'S': -1, 'T': -1, 'A': -1, 'G': -2, 'P': -2, 'D': -2, 'E': 0, 'Q': 1, 'N': 0, 'H': 0, 'R': 5, 'K': 2, 'M': -1, 'I': -3, 'L': -2, 'V': -3, 'W': -3, 'Y': -2, 'F': -3},
    'K': {'C': -3, 'S': 0, 'T': -1, 'A': -1, 'G': -2, 'P': -1, 'D': -1, 'E': 1, 'Q': 1, 'N': 0, 'H': -1, 'R': 2, 'K': 5, 'M': -1, 'I': -3, 'L': -2, 'V': -2, 'W': -3, 'Y': -2, 'F': -3},
    'M': {'C': -1, 'S': -1, 'T': -1, 'A': -1, 'G': -3, 'P': -2, 'D': -3, 'E': -2, 'Q': 0, 'N': -2, 'H': -2, 'R': -1, 'K': -1, 'M': 5, 'I': 1, 'L': 2, 'V': 1, 'W': -1, 'Y': -1, 'F': 0},
    'I': {'C': -1, 'S': -2, 'T': -1, 'A': -1, 'G': -4, 'P': -3, 'D': -3, 'E': -3, 'Q': -3, 'N': -3, 'H': -3, 'R': -3, 'K': -3, 'M': 1, 'I': 4, 'L': 2, 'V': 3, 'W': -3, 'Y': -1, 'F': 0},
    'L': {'C': -1, 'S': -2, 'T': -1, 'A': -1, 'G': -4, 'P': -3, 'D': -4, 'E': -3, 'Q': -2, 'N': -3, 'H': -3, 'R': -2, 'K': -2, 'M': 2, 'I': 2, 'L': 4, 'V': 1, 'W': -2, 'Y': -1, 'F': 0},
    'V': {'C': -1, 'S': -2, 'T': 0, 'A': 0, 'G': -3, 'P': -2, 'D': -3, 'E': -2, 'Q': -2, 'N': -3, 'H': -3, 'R': -3, 'K': -2, 'M': 1, 'I': 3, 'L': 1, 'V': 4, 'W': -3, 'Y': -1, 'F': -1},
    'W': {'C': -2, 'S': -3, 'T': -2, 'A': -3, 'G': -2, 'P': -4, 'D': -4, 'E': -3, 'Q': -2, 'N': -4, 'H': -2, 'R': -3, 'K': -3, 'M': -1, 'I': -3, 'L': -2, 'V': -3, 'W': 11, 'Y': 2, 'F': 1},
    'Y': {'C': -2, 'S': -2, 'T': -2, 'A': -2, 'G': -3, 'P': -3, 'D': -3, 'E': -2, 'Q': -1, 'N': -2, 'H': 2, 'R': -2, 'K': -2, 'M': -1, 'I': -1, 'L': -1, 'V': -1, 'W': 2, 'Y': 7, 'F': 3},
    'F': {'C': -2, 'S': -2, 'T': -2, 'A': -2, 'G': -3, 'P': -4, 'D': -3, 'E': -3, 'Q': -3, 'N': -3, 'H': -1, 'R': -3, 'K': -3, 'M': 0, 'I': 0, 'L': 0, 'V': -1, 'W': 1, 'Y': 3, 'F': 6}
}
# compute softmax of blosum62
SOFTMAX_BLOSUM62_DICT = {aa1: {aa2: np.exp(BLOSUM62_DICT[aa1][aa2]) for aa2 in AMINO_ACIDS} for aa1 in AMINO_ACIDS}
SOFTMAX_BLOSUM62_DICT = {aa1: {aa2: SOFTMAX_BLOSUM62_DICT[aa1][aa2] / sum(list(SOFTMAX_BLOSUM62_DICT[aa1].values())) for aa2 in AMINO_ACIDS} for aa1 in AMINO_ACIDS}

if __name__ == '__main__':

    os.makedirs('blosum62', exist_ok=True)

    NUM_TARGET_SAMPLES = 200

    header = "organism	mhc_class	mhc	peptide	va	ja	cdr3a	vb	jb	cdr3b	pdbid	hcnn_model"
    template_row = "human	1	A*01:01	{pep_seq}	TRAV21*01	TRAJ28*01	CAVRPGGAGPFFVVF	TRBV5-1*01	TRBJ2-7*01	CASSFNMATGQYF	5brz.pdb.human.MH1.A-01.A.C.DE	blosum62"

    for temperature in [1.0, 2.0, 3.0]:

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








