

import os, sys
import glob
import numpy as np
import pandas as pd
from scipy.special import softmax
import json
import matplotlib.pyplot as plt

from protein_holography_pytorch.utils.protein_naming import ind_to_ol_size, ol_to_ind_size

import argparse

from typing import *


'''

REMEMBER: network energies are not quite comparable scross different models and different structures

'''

def parse_pdb_for_tcrdock(path_to_tcrdock: str, pdb_dir: str, pdb: str, organism: str, mhc_class: str):

    parsed_filename = os.path.join(pdb + '__parsing_output.tsv')
    os.system(f"python {os.path.join(path_to_tcrdock, 'parse_tcr_pmhc_pdbfile.py')} \
                                --pdbfiles {os.path.join(pdb_dir, pdb+'.pdb')} \
                                --organism {organism} \
                                --mhc_class {mhc_class} \
                                --out_tsvfile {parsed_filename}")

    parsed_pdb = pd.read_csv(parsed_filename, sep='\t')
    tcrdock_info = json.loads(parsed_pdb['tcrdock_info'].values[0])
    os.remove(parsed_filename)

    organism = tcrdock_info['organism']
    mhc_class = str(tcrdock_info['mhc_class'])
    mhc = tcrdock_info['mhc_allele']
    wt_peptide = tcrdock_info['pep_seq']
    va = tcrdock_info['tcr'][0][0]
    ja = tcrdock_info['tcr'][0][1]
    cdr3a = tcrdock_info['tcr'][0][2]
    vb = tcrdock_info['tcr'][1][0]
    jb = tcrdock_info['tcr'][1][1]
    cdr3b = tcrdock_info['tcr'][1][2]

    return organism, mhc_class, mhc, wt_peptide, va, ja, cdr3a, vb, jb, cdr3b



def prepare_target_template(path_to_tcrdock: str, output_dir: str, pdb_dir: str, pdbs: List[str], organism: str, mhc_class: str):

    ## get info for tcrdock
    wt_peptides_list = []
    last_parsed_info = None

    for pdb in pdbs:
        parsed_info = parse_pdb_for_tcrdock(path_to_tcrdock, pdb_dir, pdb, organism, mhc_class)

        if last_parsed_info is None:
            last_parsed_info = parsed_info
            wt_peptides_list.append(parsed_info[3])
        else:
            # make sure all info is the same EXCEPT FOR THE PEPTIDE
            for i in range(len(last_parsed_info)):
                if i != 3:
                    assert last_parsed_info[i] == parsed_info[i]
                else:
                    assert last_parsed_info[i] != parsed_info[i]
            wt_peptides_list.append(parsed_info[3])

    organism, mhc_class, mhc, _, va, ja, cdr3a, vb, jb, cdr3b = last_parsed_info

    ## make the input tsv file with the mutant lines
    columns = ['organism', 'mhc_class', 'mhc', 'peptide', 'va', 'ja', 'cdr3a', 'vb', 'jb', 'cdr3b', 'model', 'pdbid', 'pnlogp', 'pnE', 'num_times_pep_was_sampled']
    columns_row = '\t'.join(columns)


    lines = [columns_row]
    for pdb, pep_seq in zip(pdbs, wt_peptides_list):
        line = '\t'.join([organism, str(mhc_class), mhc, pep_seq, va, ja, cdr3a, vb, jb, cdr3b, 'WT', pdb, 'nan', 'nan', 'num_times_pep_was_sampled'])
        lines.append(line)

    tcrdock_input_template_file = os.path.join(output_dir, 'tcrdock_targets_template.tsv')

    with open(tcrdock_input_template_file, 'w') as f:
        f.write('\n'.join(lines))



def get_accepted_data(data: np.ndarray, acceptance: np.ndarray):
    """Parse the accepted history from an annealing run"""
    true_data = np.zeros(
        shape=data.shape, 
        dtype=data.dtype)
    true_data[0] = data[0]
    last_datum = data[0]
    for i, (datum, acc) in enumerate(zip(data[1:], acceptance[1:])):
        if acc:
            true_data[i+1] = datum
            last_datum = datum
        else:
            true_data[i+1] = last_datum
    return true_data



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path_to_tcrdock', type=str, required=True)
    parser.add_argument('--output_dir', type=str, required=True)
    parser.add_argument('--annealing_schedule', type=str, required=True)
    parser.add_argument('--hcnn_models_to_consider', type=str, nargs='+', default=['so3_convnet_base_ensemble', 'so3_convnet_noise=0p5'])
    parser.add_argument('--pdb_dir', type=str, required=True)
    parser.add_argument('--pdbs', type=str, nargs='+', required=True)
    parser.add_argument('--organism', type=str, required=True, choices=['human', 'mouse'])
    parser.add_argument('--mhc_class', type=str, required=True, choices=['1', '2'])
    args = parser.parse_args()


    if not os.path.exists(os.path.join(args.output_dir, 'tcrdock_targets_template.tsv')):
        prepare_target_template(args.path_to_tcrdock, args.output_dir, args.pdb_dir, args.pdbs, args.organism, args.mhc_class)

    df = pd.read_csv(os.path.join(args.output_dir, 'tcrdock_targets_template.tsv'), sep='\t')

    ## add the annealing run peptides
    sequences = []
    hcnn_models = []
    template_pdbs = []
    pnlogps = []
    pnEs = []
    for base_pdb in args.pdbs:
        for hcnn_model in args.hcnn_models_to_consider:
            filenames = glob.glob(os.path.join(args.output_dir, 'annealing_runs', f'*{base_pdb}*{hcnn_model}*__{args.annealing_schedule}__*pnE*pocket*.npy'))
            for filename in filenames:
                data = np.load(filename, allow_pickle=True)[()]
                assert len(data['sequences'])-1 == len(data['peptide pnE'])
                true_sequences = get_accepted_data(data['sequences'][:-1], data['acceptance'])
                seq = true_sequences[-1].decode('utf-8')
                sequences.append(seq)
                hcnn_models.append(hcnn_model)
                template_pdbs.append(base_pdb)

                true_pes = get_accepted_data(data['peptide pes'], data['acceptance'])
                true_pnE = get_accepted_data(data['peptide pnE'], data['acceptance'])

                if len(true_pes) != len(true_pnE):
                    print(f'Error: {filename} has different lengths for pes and pnE')
                
                pes = true_pes[-1]
                pnE = true_pnE[-1]
                pnE_check = np.sum(pes[np.arange(len(pes)), [ol_to_ind_size[aa] for aa in seq]])
                assert np.isclose(pnE, pnE_check)

                logprobas = np.log(softmax(pes.astype(np.float64), axis=1))
                peptide_logproba = np.mean(logprobas[np.arange(len(logprobas)), [ol_to_ind_size[aa] for aa in seq]])
                pnlogps.append(peptide_logproba)
                pnEs.append(pnE)

    
    sequences = np.array(sequences)
    hcnn_models = np.array(hcnn_models)
    template_pdbs = np.array(template_pdbs)
    pnlogps = np.array(pnlogps)
    pnEs = np.array(pnEs)

    # collapse sequences that are the same, keep track of multiplicity though
    collapsed_sequences = []
    collapsed_pnlogps = []
    collapsed_pnEs = []
    collapsed_template_pdbs = []
    collapsed_hcnn_models = []
    collapsed_seq_multiplicites = []
    for hcnn_model in args.hcnn_models_to_consider:
        for pdb in args.pdbs:
            for seq in set(sequences[(hcnn_models == hcnn_model) & (template_pdbs == pdb)]):
                collapsed_sequences.append(seq)
                collapsed_pnlogps.append(pnlogps[(hcnn_models == hcnn_model) & (sequences == seq)].mean())
                collapsed_pnEs.append(pnEs[(hcnn_models == hcnn_model) & (sequences == seq)].mean())
                collapsed_hcnn_models.append(hcnn_model)
                collapsed_template_pdbs.append(pdb)
                collapsed_seq_multiplicites.append((sequences == seq).sum())
    
    collapsed_sequences = np.array(collapsed_sequences)
    collapsed_pnlogps = np.array(collapsed_pnlogps)
    collapsed_pnEs = np.array(collapsed_pnEs)
    collapsed_hcnn_models = np.array(collapsed_hcnn_models)
    collapsed_template_pdbs = np.array(collapsed_template_pdbs)
    collapsed_seq_multiplicites = np.array(collapsed_seq_multiplicites)

    # save sequences to df, one per row
    # copy the previous row for each sequence, then overwrite the sequence
    template_row = df.iloc[0]

    for pep_seq, hcnn_model, template_pdb, pnlogp, pnE, pep_mult in zip(collapsed_sequences, collapsed_hcnn_models, collapsed_template_pdbs, collapsed_pnlogps, collapsed_pnEs, collapsed_seq_multiplicites):
        new_row = template_row.copy()
        new_row['peptide'] = pep_seq
        new_row['pdbid'] = template_pdb
        new_row['model'] = hcnn_model
        new_row['pnlogp'] = pnlogp
        new_row['pnE'] = pnE
        new_row['num_times_pep_was_sampled'] = pep_mult
        df = pd.concat([df, new_row.to_frame().T], ignore_index=True)
    
    df.to_csv(f'tcrdock_targets_with_sampled_peptides__{args.annealing_schedule}__{"__".join(args.hcnn_models_to_consider)}.tsv', sep='\t', index=False)








