
import os
import numpy as np
import pandas as pd
import gzip, pickle

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
os.makedirs('plots', exist_ok=True)

from compare_pwms import normalize_pwm, entropy_of_pwm, kl_divergence_of_pwms


import sys
sys.path.append('../mutation_effects/src')
from compute_rmsd import extract_chain_residue_info, build_chain_sequences, align_sequences_middle_gap_penalty, count_initial_gap, align_chains, get_ca_atoms, compute_rmsd



# every PWM here follows this numbering
# from hermes.utils.protein_naming import ind_to_ol_size, ol_to_ind_size
ind_to_ol_size = {0: 'G', 1: 'A', 2: 'C', 3: 'S', 4: 'P', 5: 'T', 6: 'V', 7: 'D', 8: 'I', 9: 'L', 10: 'N', 11: 'M', 12: 'Q', 13: 'K', 14: 'E', 15: 'H', 16: 'F', 17: 'R', 18: 'Y', 19: 'W'}
ol_to_ind_size = {ind_to_ol_size[key]: key for key in ind_to_ol_size}

def get_individual_pwm_path(df_row):
    return f"pwm_csv_files/mhc_crystal_hcnn_fixed_structure/hermes_py_000/{df_row['pdbid']}__{df_row['wt_peptide']}__{df_row['mhc_allele']}.csv"

def get_pdbfile(df_row):
    return f"pdbs/pmhc/{df_row['pdbid']}.pdb"

def get_mhc_chain(seq_dict, pep_chain):
    # longest chain that is not peptide chain is most likely the MHC chain (if there is something else other than MHC, like a ligand, it will most likely have a shorter chain)
    chains = []
    lengths = []
    for chain, seq in seq_dict.items():
        if chain != pep_chain:
            chains.append(chain)
            lengths.append(len(seq))
    return chains[np.argmax(lengths)]

def align_on_mhc_and_compute_rmsds(pdbfile_1, pep_chain_1, pdbfile_2, pep_chain_2):

    # Extract chain information
    chain_dict_1 = extract_chain_residue_info(pdbfile_1)
    seq_dict_1 = build_chain_sequences(chain_dict_1)
    mhc_chain_1 = get_mhc_chain(seq_dict_1, pep_chain_1)

    chain_dict_2 = extract_chain_residue_info(pdbfile_2)
    seq_dict_2 = build_chain_sequences(chain_dict_2)
    mhc_chain_2 = get_mhc_chain(seq_dict_2, pep_chain_2)

    # Align MHC chains
    aligned_seq_mhc_1, aligned_seq_mhc_2, score = align_sequences_middle_gap_penalty(seq_dict_1[mhc_chain_1], seq_dict_2[mhc_chain_2], mode='global')
    # print(f"Aligned MHC Sequences:\n{aligned_seq_mhc_1}\n{aligned_seq_mhc_2}")
    initial_gap_1 = count_initial_gap(aligned_seq_mhc_1)
    initial_gap_2 = count_initial_gap(aligned_seq_mhc_2)
    common_indices = np.array([i for i, (a, b) in enumerate(zip(aligned_seq_mhc_1, aligned_seq_mhc_2)) if a != '-' and b != '-'])
    resnums_mhc_1 = np.array(chain_dict_1[mhc_chain_1]['resnums'])[common_indices - initial_gap_1]
    resnums_mhc_2 = np.array(chain_dict_2[mhc_chain_2]['resnums'])[common_indices - initial_gap_2]
    _, structure1, structure2 = align_chains(pdbfile_1, mhc_chain_1, resnums_mhc_1,
                                                pdbfile_2, mhc_chain_2, resnums_mhc_2)
    mhc_atoms_1 = get_ca_atoms(structure1, mhc_chain_1, resnums_mhc_1)
    mhc_atoms_2 = get_ca_atoms(structure2, mhc_chain_2, resnums_mhc_2)
    mhc_rmsd = compute_rmsd(mhc_atoms_1, mhc_atoms_2)
    print(f'RMSD of MHCs aligning MHCs: {mhc_rmsd:.2f}')

    # # compute RMSD for peptide chain aligned on the MHC
    # aligned_seq_pep_1, aligned_seq_pep_2, score = align_sequences_middle_gap_penalty(seq_dict_1[pep_chain_1], seq_dict_2[pep_chain_2], mode='global')
    # print(f"Aligned Peptide Sequences:\n{aligned_seq_pep_1}\n{aligned_seq_pep_2}")
    # initial_gap_1 = count_initial_gap(aligned_seq_pep_1)
    # initial_gap_2 = count_initial_gap(aligned_seq_pep_2)
    # common_indices = np.array([i for i, (a, b) in enumerate(zip(aligned_seq_pep_1, aligned_seq_pep_2)) if a != '-' and b != '-'])
    # resnums_pep_1 = np.array(chain_dict_1[pep_chain_1]['resnums'])[common_indices - initial_gap_1]
    # resnums_pep_2 = np.array(chain_dict_2[pep_chain_2]['resnums'])[common_indices - initial_gap_2]
    resnums_pep_1 = np.array(chain_dict_1[pep_chain_1]['resnums'])
    resnums_pep_2 = np.array(chain_dict_2[pep_chain_2]['resnums'])
    pep_atoms_1 = get_ca_atoms(structure1, pep_chain_1, resnums_pep_1)
    pep_atoms_2 = get_ca_atoms(structure2, pep_chain_2, resnums_pep_2)
    pep_rmsd = compute_rmsd(pep_atoms_1, pep_atoms_2)
    print(f'RMSD of Peptides aligning MHCs: {pep_rmsd:.2f}')
    print()

    return pep_rmsd


def get_pwms(final_pwms, df, allele, length):

    df_rows = df.loc[np.logical_and(df['mhc_allele_final'] == allele, df['peptide_length'] == peptide_length)]

    assert len(df_rows) == 2, (allele, length, len(df_rows), final_pwms[allele][length]['num_structures'])

    df_row_1 = df_rows.iloc[0]
    df_row_2 = df_rows.iloc[1]

    pwm_1 = pd.read_csv(get_individual_pwm_path(df_row_1)).values
    pwm_2 = pd.read_csv(get_individual_pwm_path(df_row_2)).values

    pwm_combined = final_pwms[allele][length]['hermes_py_000_pwm']
    pwm_true = final_pwms[allele][length]['mhc_motif_pwm']

    return pwm_1, pwm_2, pwm_combined, pwm_true


def make_distribution_plot(metrics, ylabel, outpath):

    conf_to_color = {
        'same_conf': 'blue',
        'diff_conf': 'orange'
    }

    fontsize = 15

    jitter_strength = 0.08
    mean_line_thickness = 1.2

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(4, 4))

    # Plotting loop
    conf_list = list(metrics.keys())
    for i, conf in conf_list:

        data = metrics[conf]
        color = conf_to_color[conf]

        x_center = i + 1
        x_jitter = x_center + np.random.uniform(-jitter_strength, jitter_strength, size=len(data))
        
        # Scatter plot with jitter
        ax.scatter(x_jitter, data, color=color, alpha=0.7, label=conf)
        
        # Plot mean line
        mean_val = np.mean(data)
        ax.hlines(mean_val, x_center - 0.2, x_center + 0.2, colors=color, linewidth=mean_line_thickness)

    # Customize x-axis
    ax.set_xticks([1, 2])
    ax.set_xticklabels(conf_list, fontsize=fontsize)

    ax.tick_params(axis='y', labelsize=fontsize-2)

    ax.set_ylabel(ylabel, fontsize=fontsize)

    # Add legend
    ax.legend()

    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


df = pd.read_csv('pmhc_class_1_crystal_structures.csv')

with gzip.open('pmhc_class_1_crystal_structures_pwms.pkl.gz', 'rb') as f:
    final_pwms = pickle.load(f)

allele_and_length_list = []
pep_rmsd_list = []

for allele in final_pwms:
    for peptide_length in final_pwms[allele]:

        allele_pwms = final_pwms[allele][peptide_length] 

        if 'num_structures' not in allele_pwms:
            continue

        if allele_pwms['num_structures'] == 2:

            df_rows = df.loc[np.logical_and(df['mhc_allele_final'] == allele, df['peptide_length'] == peptide_length)]

            if len(df_rows) == 2:

                print((allele, peptide_length, len(df_rows)))

                df_row_1 = df_rows.iloc[0]
                df_row_2 = df_rows.iloc[1]

                pdbfile_1 = get_pdbfile(df_row_1)
                pdbfile_2 = get_pdbfile(df_row_2)

                pep_rmsd = align_on_mhc_and_compute_rmsds(pdbfile_1, df_row_1['peptide_chain'], pdbfile_2, df_row_2['peptide_chain'])

                allele_and_length_list.append((allele, peptide_length))
                pep_rmsd_list.append(pep_rmsd)

allele_and_length_list = np.array(allele_and_length_list)
pep_rmsd_list = np.array(pep_rmsd_list)

pep_rmsd_cutoff = 1.5

allele_and_length_list__same_conf = allele_and_length_list[pep_rmsd_list <= pep_rmsd_cutoff]
allele_and_length_list__diff_conf = allele_and_length_list[pep_rmsd_list > pep_rmsd_cutoff]
print(len(allele_and_length_list))
print(len(allele_and_length_list__same_conf))
print(len(allele_and_length_list__diff_conf))


all_metrics = {'kld_percent_improvement': {'same_conf': [], 'diff_conf': []},
           'symmetric_kld_between_the_two_pwms': {'same_conf': [], 'diff_conf': []}}

for curr__allele_and_length_list, conf in [(allele_and_length_list__same_conf, 'same_conf'), (allele_and_length_list__diff_conf, 'diff_conf')]:
                                           
    for (allele, length) in curr__allele_and_length_list:

        pwm_1, pwm_2, pwm_combined, pwm_true = get_pwms(final_pwms, df, allele, int(length))

        kld_1 = kl_divergence_of_pwms(pwm_true, pwm_1) / length
        kld_2 = kl_divergence_of_pwms(pwm_true, pwm_2) / length
        kld_combined = kl_divergence_of_pwms(pwm_true, pwm_combined) / length
        kld_percent_improvement = ((kld_combined - np.mean([kld_1, kld_2])) / np.mean([kld_1, kld_2])) * 100

        symmetric_kld_between_the_two_pwms = ((kl_divergence_of_pwms(pwm_1, pwm_2) + kl_divergence_of_pwms(pwm_2, pwm_1)) / 2) / length

        all_metrics['kld_percent_improvement'][conf].append(kld_percent_improvement)
        all_metrics['symmetric_kld_between_the_two_pwms'][conf].append(symmetric_kld_between_the_two_pwms)


make_distribution_plot(all_metrics['kld_percent_improvement'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_kld_percent_improvement.png')
make_distribution_plot(all_metrics['symmetric_kld_between_the_two_pwms'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_symmetric_kld_between_the_two_pwms.png')






                         





