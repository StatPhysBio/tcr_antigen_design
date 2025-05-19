
import os
import gzip, pickle

import numpy as np
import pandas as pd
from scipy import stats

from tqdm import tqdm

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
os.makedirs('plots', exist_ok=True)

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

from compare_pwms import normalize_pwm, entropy_of_pwm, kl_divergence_of_pwms


import sys
sys.path.append('../mutation_effects/src')
from compute_rmsd import extract_chain_residue_info, build_chain_sequences, align_sequences_middle_gap_penalty, count_initial_gap, align_chains, get_ca_atoms, compute_rmsd



# every PWM here follows this numbering
# from hermes.utils.protein_naming import ind_to_ol_size, ol_to_ind_size
ind_to_ol_size = {0: 'G', 1: 'A', 2: 'C', 3: 'S', 4: 'P', 5: 'T', 6: 'V', 7: 'D', 8: 'I', 9: 'L', 10: 'N', 11: 'M', 12: 'Q', 13: 'K', 14: 'E', 15: 'H', 16: 'F', 17: 'R', 18: 'Y', 19: 'W'}
ol_to_ind_size = {ind_to_ol_size[key]: key for key in ind_to_ol_size}

def make_pretty_allele(allele):
    if '-' in allele: # mouse
        pretty_allele = allele
    else:
        pretty_allele = allele[0] + '*' + allele[1:3] + ':' + allele[3:]
    return pretty_allele

def hamming_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
    return sum(a != b for a, b in zip(seq1, seq2))

def get_individual_pwm_path(df_row):
    return f"pwm_csv_files/mhc_crystal_hcnn_fixed_structure/hermes_py_000/{df_row['pdbid'].values[0]}__{df_row['wt_peptide'].values[0]}__{df_row['mhc_allele'].values[0]}.csv"

def get_pdbfile(df_row):
    return f"pdbs/pmhc/{df_row['pdbid'].values[0]}.pdb"

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
    # print(f'RMSD of MHCs aligning MHCs: {mhc_rmsd:.2f}')

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
    # print(f'RMSD of Peptides aligning MHCs: {pep_rmsd:.2f}')
    # print()

    return mhc_rmsd, pep_rmsd


def get_pwms(final_pwms, df, allele, peptide_length, pdb_pair):

    pdbs = final_pwms[allele][peptide_length]['pdbs'][pdb_pair[0]], final_pwms[allele][peptide_length]['pdbs'][pdb_pair[1]]

    df_row_1 = df.loc[np.logical_and(np.logical_and(df['mhc_allele_final'] == allele, df['peptide_length'] == peptide_length), df['pdbid'] == pdbs[0])]
    df_row_2 = df.loc[np.logical_and(np.logical_and(df['mhc_allele_final'] == allele, df['peptide_length'] == peptide_length), df['pdbid'] == pdbs[1])]

    pwm_1 = pd.read_csv(get_individual_pwm_path(df_row_1), index_col=0).values
    pwm_2 = pd.read_csv(get_individual_pwm_path(df_row_2), index_col=0).values

    # pwm_combined = final_pwms[allele][peptide_length]['hermes_py_000_pwm']
    pwm_combined = normalize_pwm((pwm_1 + pwm_2) / 2)
    pwm_true = final_pwms[allele][peptide_length]['mhc_motif_pwm']

    pep_1 = df_row_1['wt_peptide'].values[0]
    pep_2 = df_row_2['wt_peptide'].values[0]

    return pwm_1, pwm_2, pwm_combined, pwm_true, pep_1, pep_2, pdbs[0], pdbs[1]


def hamming_distance_to_marker_size(hdist):
    return list(np.log(((np.array(hdist) + 1)*1.1))*200)



df = pd.read_csv('pmhc_class_1_crystal_structures.csv', index_col=0)

with gzip.open('pmhc_class_1_crystal_structures_pwms.pkl.gz', 'rb') as f:
    final_pwms = pickle.load(f)

from itertools import combinations

np.random.seed(42)

# max_num_pairs_per_mhc = 500

# allele_and_length_list = []
# pep_rmsd_list = []
# kld_over_l_list = []
# kld_single_pwm_list = []
# num_structures_list = []

# for allele in tqdm(final_pwms):
#     for peptide_length in final_pwms[allele]:

#         allele_pwms = final_pwms[allele][peptide_length]

#         if 'num_structures' not in allele_pwms or 'hermes_py_000_pwm' not in allele_pwms or 'mhc_motif_pwm' not in allele_pwms:
#             continue

#         if allele_pwms['num_structures'] >= 2:

#             kld_over_l = kl_divergence_of_pwms(allele_pwms['hermes_py_000_pwm'], allele_pwms['mhc_motif_pwm']) / peptide_length

#             combos = list(combinations(np.arange(allele_pwms['num_structures']), 2))
#             np.random.shuffle(combos)

#             pep_rmsd_temp_list = []
#             kld_single_pwm_temp_list = []
#             idxs_done = set()
            
#             for combo in tqdm(combos[:max_num_pairs_per_mhc]):

#                 idx1, idx2 = combo

#                 pdbs = allele_pwms['pdbs'][idx1], allele_pwms['pdbs'][idx2]

#                 df_row_1 = df.loc[np.logical_and(np.logical_and(df['mhc_allele_final'] == allele, df['peptide_length'] == peptide_length), df['pdbid'] == pdbs[0])]
#                 df_row_2 = df.loc[np.logical_and(np.logical_and(df['mhc_allele_final'] == allele, df['peptide_length'] == peptide_length), df['pdbid'] == pdbs[1])]

#                 assert len(df_row_1) == 1, df_row_1
#                 assert len(df_row_2) == 1, len(df_row_2)

#                 pdbfile_1 = get_pdbfile(df_row_1)
#                 pdbfile_2 = get_pdbfile(df_row_2)

#                 pwm_1 = pd.read_csv(get_individual_pwm_path(df_row_1), index_col=0).values
#                 pwm_2 = pd.read_csv(get_individual_pwm_path(df_row_2), index_col=0).values

#                 for idx, pwm in zip(combo, [pwm_1, pwm_2]):
#                     if idx not in idxs_done:
#                         kld_single_pwm_temp_list.append(kl_divergence_of_pwms(pwm, allele_pwms['mhc_motif_pwm']) / peptide_length)
#                         idxs_done.add(idx)
                
#                 try:
#                     mhc_rmsd, pep_rmsd = align_on_mhc_and_compute_rmsds(pdbfile_1, df_row_1['peptide_chain'].values[0], pdbfile_2, df_row_2['peptide_chain'].values[0])
#                 except Exception as e:
#                     print('Error in alignment step: ', e)
#                     continue

#                 # print('success')

#                 if mhc_rmsd < 1.5: # otherwise we shouldn't trust the alignment
#                     if pep_rmsd < 30: # I think these are buggy, just a handful of them
#                         pep_rmsd_temp_list.append(pep_rmsd)
            
#             if len(pep_rmsd_temp_list) > 0:
#                 allele_and_length_list.append((allele, peptide_length))
#                 pep_rmsd_list.append(np.mean(pep_rmsd_temp_list))
#                 kld_over_l_list.append(kld_over_l)
#                 kld_single_pwm_list.append(np.mean(kld_single_pwm_temp_list))
#                 num_structures_list.append(allele_pwms['num_structures'])

# results = {'pep_rmsd_list': pep_rmsd_list,
#                  'allele_and_length_list': allele_and_length_list,
#                  'kld_over_l_list': kld_over_l_list,
#                  'kld_single_pwm_list': kld_single_pwm_list,
#                  'num_structures_list': num_structures_list}

# with gzip.open('pmhc_rmsd_and_kld.pkl.gz', 'wb') as f:
#     pickle.dump(results, f)



with gzip.open('pmhc_rmsd_and_kld.pkl.gz', 'rb') as f:
    results = pickle.load(f)

peptide_length_to_color = {
    8: 'yellow',
    9: 'tab:orange',
    10: 'tab:red',
    11: 'tab:purple'
}

def allele_to_marker(allele):
    # circle for human, square for mouse
    if '-' in allele:
        return 'd'
    else:
        return 'o'

def num_structures_to_marker_size(num_structures_list):
    return list(np.log(((np.array(num_structures_list) + 1)*1.05))*100)

def rmsd_to_marker_size(rmsd_list):
    return list(np.log(((np.array(rmsd_list) + 1)*1.07))*150)

xs = results['num_structures_list']
data = np.array(results['kld_over_l_list']) # / np.array(results['kld_single_pwm_list'])
colors = [peptide_length_to_color[pep_length] for _, pep_length in results['allele_and_length_list']]
sizes = rmsd_to_marker_size(results['pep_rmsd_list'])
markers = [allele_to_marker(allele) for allele, _ in results['allele_and_length_list']]


fontsize = 20

# Create figure and axis
fig, ax = plt.subplots(figsize=(4.5, 4))

# Scatter plot with jitter
for x, d, c, s, m in zip(xs, data, colors, sizes, markers):
    ax.scatter(x, d, facecolor=c, alpha=0.7, edgecolor='black', marker=m, s=140)


sr, sr_pval = stats.spearmanr(xs, data)

ax.text(0.95, 0.90, rf'$\rho$: {sr:.2f}', ha='right', fontsize=fontsize, transform=ax.transAxes)
# ax.text(0.95, 0.80, f'p-val: {sr_pval:.6f}', ha='right', fontsize=fontsize-1, transform=ax.transAxes)


ax.tick_params(axis='both', labelsize=fontsize)

# ax.set_xlabel('mean peptide RMSD', fontsize=fontsize)
ax.set_xscale('log')
ax.set_xlabel('number of pMHC structures', fontsize=fontsize)
ax.set_ylabel('KL-Divergence / L\nto MHC motif PWM', fontsize=fontsize)

# # Add legend
# ax.legend()

plt.tight_layout()
plt.savefig('pmhc_pwm_peptide_rmsd_vs_kl_divergence.png', bbox_inches='tight')
plt.savefig('pmhc_pwm_peptide_rmsd_vs_kl_divergence.pdf', bbox_inches='tight')
plt.close()


legend_handles = [
    Line2D([0], [0], marker='o', color='black', linestyle='None', markersize=8, label='human'),
    Line2D([0], [0], marker='d', color='black', linestyle='None', markersize=8, label='mouse')
]
fig, ax = plt.subplots(figsize=(3, 3))
ax.axis('off')
ax.legend(handles=legend_handles, fontsize=13, title_fontsize=14)
plt.savefig('legend__organism_markers.png')
plt.savefig('legend__organism_markers.pdf')
plt.close()


legend_handles = [
    Patch(facecolor=color, edgecolor='black', label=length)
    for length, color in peptide_length_to_color.items()
]
fig, ax = plt.subplots(figsize=(3, 3))
ax.axis('off')
ax.legend(handles=legend_handles, title='peptide length', fontsize=13, title_fontsize=14, ncol=2)
plt.savefig('legend__peptide_length_colors.png')
plt.savefig('legend__peptide_length_colors.pdf')
plt.close()



