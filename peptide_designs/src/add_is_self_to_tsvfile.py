import os
import numpy as np
import pandas as pd
import time
from tqdm import tqdm
import argparse
from Bio import SeqIO

'''

1. Adds column 'is_self_peptide' to the input .tsv file
2. Adds column 'self_protein_name' to the input .tsv file
3. Saves the new .tsv file with the same name, unless output_tsvfile is specified

'''

PATH_TO_FASTA_OF_SELF_PROTEINOME = '/gscratch/spe/gvisan01/self_peptides_lookup/UP000005640_9606.fasta'

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_tsvfile', type=str, required=True)
    parser.add_argument('-o', '--output_tsvfile', type=str, default=None)
    parser.add_argument('-p', '--peptides_column', type=str, default='peptide')
    args = parser.parse_args()

    if args.output_tsvfile is None:
        output_tsvfile = args.input_tsvfile
    else:
        output_tsvfile = args.output_tsvfile

    df = pd.read_csv(args.input_tsvfile, sep='\t')

    peptides = df[args.peptides_column].values



    fasta_sequences = SeqIO.parse(open(PATH_TO_FASTA_OF_SELF_PROTEINOME), 'fasta')

    def is_in_self_brute_force(peptide, fasta_sequences):
        for fasta in fasta_sequences:
            name, sequence = fasta.id, str(fasta.seq)
            if peptide in sequence:
                print(name)
                return 1, name
        return 0, None

    is_self_peptide, self_protein_name = [], []
    for peptide in tqdm(peptides):
        is_self, name = is_in_self_brute_force(peptide, fasta_sequences)
        is_self_peptide.append(is_self)
        self_protein_name.append(name)

    df['is_self_peptide'] = is_self_peptide
    df['self_protein_name'] = self_protein_name

    df.to_csv(output_tsvfile, sep='\t', index=False)
    