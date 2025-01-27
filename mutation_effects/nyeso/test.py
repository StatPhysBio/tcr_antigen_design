
import numpy as np
import pandas as pd

df = pd.read_csv('nyeso_peptide_kd_closest.csv', index_col=False)

peptides = df['sequence'].values

print(df)

print(peptides)

def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    return sum([seq1[i] != seq2[i] for i in range(len(seq1))])

all_hamming_distances = np.array([[hamming_distance(seq1, seq2) for seq1 in peptides] for seq2 in peptides])

print(np.max(all_hamming_distances))
