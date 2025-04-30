
import os
import numpy as np
import pandas as pd
import argparse


ESM_PATH = '/gscratch/spe/gvisan01/esm/esm'



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--pdbdir', type=str, required=True)
    parser.add_argument('--csv_filepath', type=str, required=True)
    parser.add_argument('--pdb_column', type=str, required=True)
    parser.add_argument('--chain_column', type=str, required=True)
    parser.add_argument('--peptide_column', type=str, required=True)
    parser.add_argument('--output_csv_filepath', type=str, required=True)
    args = parser.parse_args()

    output_directory = os.path.dirname(args.output_csv_filepath)
    os.makedirs(output_directory, exist_ok=True)

    # make temporary directory with fasta files and scores
    temp_fasta_dir = os.path.join(output_directory, 'temp_fasta')
    os.makedirs(temp_fasta_dir, exist_ok=True)
    temp_output_dir = os.path.join(output_directory, 'temp_output')
    os.makedirs(temp_output_dir, exist_ok=True)

    # load the csv file
    df = pd.read_csv(args.csv_filepath)

    # initialize peptide->score dictionary
    peptide_to_score = {}

    # group the df by (pdb, chain) tuples
    grouped_df = df.groupby([args.pdb_column, args.chain_column])

    for df_group in grouped_df:
        # get the pdb and chain
        pdb, chain = df_group[0]
        # get the peptide sequences
        peptides = df_group[1][args.peptide_column].values

        # make fasta file, call it with pdb and chain
        fasta_file_path = os.path.join(temp_fasta_dir, f'{pdb}_{chain}.fasta')
        with open(fasta_file_path, 'w') as fasta_file:
            for i, peptide in enumerate(peptides):
                fasta_file.write(f'>{peptide}\n{peptide}\n')
        
        # score with ESM-IF1
        esmif_output_file_path = os.path.join(temp_output_dir, f'{pdb}_{chain}_scores.csv')

        command = f'python {ESM_PATH}/examples/inverse_folding/score_log_likelihoods.py {os.path.join(args.pdbdir, pdb+".pdb")} {fasta_file_path} --chain {chain} --outpath {esmif_output_file_path} --multichain-backbone'

        print(command)

        os.system(command)
    
        # write the scores in a peptide->score dictionary
        scores_df = pd.read_csv(esmif_output_file_path)
        for i, row in scores_df.iterrows():
            peptide = row['seqid']
            score = row['log_likelihood']
            peptide_to_score[peptide] = score
        
    # write the scores in the original dataframe
    df['esmif_log_proba'] = df.apply(lambda row: peptide_to_score.get(row[args.peptide_column], np.nan), axis=1)

    # save the dataframe
    df.to_csv(args.output_csv_filepath, index=False)





