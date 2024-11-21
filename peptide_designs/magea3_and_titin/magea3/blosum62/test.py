
import os
import numpy as np
import pandas as pd



if __name__ == '__main__':

    df = pd.read_csv('sample_peptides_from_blosum62__temperature=1.0.tsv', sep='\t', index_col=False)
    # df = pd.read_csv('../mhc_pwm/mhc_motif_peptides.tsv', sep='\t')

    print(df.head())
