
import os
import gzip, pickle

import numpy as np
import pandas as pd
from scipy import stats

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

def make_distribution_plot(metrics, ylabel, outpath, conf_to_hdist, conf_to_color=None):

    conf_to_pretty_name = {
        'same_conf': 'Same\nConf.',
        'diff_conf': 'Different\nConf.'        
    }

    conf_to_marker = {
        'same_conf': 'o',
        'diff_conf': 'o'
    }

    np.random.seed(42)

    fontsize = 17

    jitter_strength = 0.15
    mean_line_thickness = 1.5

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(4, 4))

    # Plotting loop
    conf_list = list(metrics.keys())
    for i, conf in enumerate(conf_list):

        data = metrics[conf]
        colors = conf_to_color[conf]
        sizes = hamming_distance_to_marker_size(conf_to_hdist[conf])
        marker = conf_to_marker[conf]

        x_center = i + 1
        x_jitter = x_center + np.random.uniform(-jitter_strength, jitter_strength, size=len(data))
        
        # Scatter plot with jitter
        for x, d, c, s in zip(x_jitter, data, colors, sizes):
            ax.scatter(x, d, facecolor=c, alpha=0.7, s=s, edgecolor='black', marker=marker)
        
        # Plot mean line
        mean_val = np.mean(data)
        ax.hlines(mean_val, x_center - 0.2, x_center + 0.2, colors='black', linewidth=mean_line_thickness)
    
    t_stat, p_value = stats.ttest_ind(metrics[conf_list[0]], metrics[conf_list[1]], equal_var=False)
    
    ax.text(0.50, 0.98, f'ttest p-val: {p_value:.2f}', ha='center', va='top', fontsize=fontsize-2, transform=ax.transAxes)
    # ax.set_title( f'p-val: {p_value:.2f}', fontsize=fontsize)

    # Customize x-axis
    ax.set_xticks([1, 2])
    ax.set_xticklabels([conf_to_pretty_name[conf] for conf in conf_list], fontsize=fontsize, rotation=70)

    ax.tick_params(axis='y', labelsize=fontsize)

    ax.set_ylabel(ylabel, fontsize=fontsize)

    # # Add legend
    # ax.legend()

    plt.tight_layout()
    plt.savefig(outpath, bbox_inches='tight')
    plt.savefig(outpath.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()


def make_correlation_plot(metrics, ylabel, outpath, conf_to_x, conf_to_size_base, xlabel, conf_to_color=None):

    conf_to_pretty_name = {
        'same_conf': 'Same\nConf.',
        'diff_conf': 'Different\nConf.'        
    }

    conf_to_marker = {
        'same_conf': 'o',
        'diff_conf': 'o'
    }

    fontsize = 17

    # Create figure and axis
    fig, ax = plt.subplots(figsize=(4, 3.5))


    all_xs = []
    all_data = []    
    conf_list = list(metrics.keys())
    for i, conf in enumerate(conf_list):

        data = metrics[conf]
        colors = conf_to_color[conf]
        xs = conf_to_x[conf]
        sizes = hamming_distance_to_marker_size(conf_to_size_base[conf])
        marker = conf_to_marker[conf]

        # Scatter plot with jitter
        for x, d, c, s in zip(xs, data, colors, sizes):
            ax.scatter(x, d, facecolor=c, alpha=0.7, s=s, edgecolor='black', marker=marker)

        all_xs.append(xs)
        all_data.append(data)

    all_xs = np.hstack(all_xs)
    all_data = np.hstack(all_data)
    pr, pr_pval = stats.pearsonr(all_xs, all_data)

    ax.set_title(f'R: {pr:.2f}, p-val: {pr_pval:.2e}', fontsize=fontsize-1)
    # ax.text(0.95, 0.14, f'R: {pr:.2f}', ha='right', fontsize=fontsize-2, transform=ax.transAxes)
    # ax.text(0.95, 0.05, f'p-val: {pr_pval:.2e}', ha='right', fontsize=fontsize-2, transform=ax.transAxes)


    ax.tick_params(axis='both', labelsize=fontsize)

    ax.set_xlabel(xlabel, fontsize=fontsize)
    ax.set_ylabel(ylabel, fontsize=fontsize)

    # # Add legend
    # ax.legend()

    plt.tight_layout()
    plt.savefig(outpath, bbox_inches='tight')
    plt.savefig(outpath.replace('.png', '.pdf'), bbox_inches='tight')
    plt.close()



df = pd.read_csv('pmhc_class_1_crystal_structures.csv', index_col=0)

with gzip.open('pmhc_class_1_crystal_structures_pwms.pkl.gz', 'rb') as f:
    final_pwms = pickle.load(f)

from itertools import combinations

np.random.seed(42)

max_num_pairs_per_mhc = 20

allele_and_length_list = []
pep_rmsd_list = []
pdb_indices = []

for allele in final_pwms:
    for peptide_length in final_pwms[allele]:

        allele_pwms = final_pwms[allele][peptide_length] 

        if 'num_structures' not in allele_pwms:
            continue

        if allele_pwms['num_structures'] >= 2:

            combos = list(combinations(np.arange(allele_pwms['num_structures']), 2))
            np.random.shuffle(combos)
            
            for combo in combos[:max_num_pairs_per_mhc]:

                # print(combo)
                # print(allele, peptide_length)

                pdbs = allele_pwms['pdbs'][combo[0]], allele_pwms['pdbs'][combo[1]]

                df_row_1 = df.loc[np.logical_and(np.logical_and(df['mhc_allele_final'] == allele, df['peptide_length'] == peptide_length), df['pdbid'] == pdbs[0])]
                df_row_2 = df.loc[np.logical_and(np.logical_and(df['mhc_allele_final'] == allele, df['peptide_length'] == peptide_length), df['pdbid'] == pdbs[1])]

                assert len(df_row_1) == 1, df_row_1
                assert len(df_row_2) == 1, len(df_row_2)

                pdbfile_1 = get_pdbfile(df_row_1)
                pdbfile_2 = get_pdbfile(df_row_2)

                try:
                    mhc_rmsd, pep_rmsd = align_on_mhc_and_compute_rmsds(pdbfile_1, df_row_1['peptide_chain'].values[0], pdbfile_2, df_row_2['peptide_chain'].values[0])
                except Exception as e:
                    print('Error in alignment step: ', e)
                    continue

                if mhc_rmsd < 1.5: # otherwise we shouldn't trust the alignment
                    if pep_rmsd < 30: # I think these are buggy, just a handful of them
                        allele_and_length_list.append((allele, peptide_length))
                        pep_rmsd_list.append(pep_rmsd)
                        pdb_indices.append(combo)


allele_and_length_list = np.array(allele_and_length_list)
pep_rmsd_list = np.array(pep_rmsd_list)
pdb_indices = np.array(pdb_indices)

pep_rmsd_cutoff = 1.5

allele_and_length_list__same_conf = allele_and_length_list[pep_rmsd_list <= pep_rmsd_cutoff]
allele_and_length_list__diff_conf = allele_and_length_list[pep_rmsd_list > pep_rmsd_cutoff]

pep_rmsd_list__same_conf = pep_rmsd_list[pep_rmsd_list <= pep_rmsd_cutoff]
pep_rmsd_list__diff_conf = pep_rmsd_list[pep_rmsd_list > pep_rmsd_cutoff]

pdb_indices__same_conf = pdb_indices[pep_rmsd_list <= pep_rmsd_cutoff]
pdb_indices__diff_conf = pdb_indices[pep_rmsd_list > pep_rmsd_cutoff]

print(len(allele_and_length_list))
print(len(allele_and_length_list__same_conf))
print(len(allele_and_length_list__diff_conf))


all_metrics = {'kld_absolute_improvement': {'same_conf': [], 'diff_conf': []},
                'kld_percent_improvement': {'same_conf': [], 'diff_conf': []},
                'symmetric_kld_between_the_two_pwms': {'same_conf': [], 'diff_conf': []},
                'hamming_distance_between_the_two_peptides': {'same_conf': [], 'diff_conf': []},
                'peptide_length': {'same_conf': [], 'diff_conf': []},
                'marker': {'same_conf': [], 'diff_conf': []},
                'color': {'same_conf': [], 'diff_conf': []},
                'pep_rmsd': {'same_conf': [], 'diff_conf': []}}

peptide_length_to_marker = {
    8: 'o',
    9: 's',
    10: 'd',
    11: 'x'
}

peptide_length_to_color = {
    8: 'yellow',
    9: 'tab:orange',
    10: 'tab:red',
    11: 'tab:purple'
}

from matplotlib.patches import Patch
from matplotlib.lines import Line2D

legend_handles = [
    Patch(facecolor=color, edgecolor='black', label=length)
    for length, color in peptide_length_to_color.items()
]
fig, ax = plt.subplots(figsize=(3, 3))
ax.axis('off')
ax.legend(handles=legend_handles, title='peptide length', fontsize=13, title_fontsize=14)
plt.savefig('plots/legend__peptide_length_colors.png')
plt.savefig('plots/legend__peptide_length_colors.pdf')
plt.close()

conf_to_marker = {
    'same conf.': 'o',
    'different conf.': 's'
}

legend_handles = [
    Line2D([0], [0], marker=marker, color='black', linestyle='None',
           markersize=8, label=label)
    for label, marker in conf_to_marker.items()
]
fig, ax = plt.subplots(figsize=(3, 3))
ax.axis('off')
ax.legend(handles=legend_handles, fontsize=13, title_fontsize=14)
plt.savefig('plots/legend__conf_markers.png')
plt.savefig('plots/legend__conf_markers.pdf')
plt.close()

hdist_to_show = [0, 0.25, 0.5, 0.75, 1]
sizes = hamming_distance_to_marker_size(hdist_to_show)

legend_handles = [
    Line2D([0], [0], marker='o', color='black', linestyle='None',
           markersize=size**0.5,  # convert area to diameter
           label=f'{hdist:.2f}')
    for hdist, size in zip(hdist_to_show, sizes)
]

fig, ax = plt.subplots(figsize=(3, 3))
ax.axis('off')
ax.legend(handles=legend_handles, title='hamming\ndistance (%)', fontsize=13, title_fontsize=14)
plt.savefig('plots/legend__hdist_size.png')
plt.savefig('plots/legend__hdist_size.pdf')
plt.close()

table_for_paper = {'MHC allele': [], 'peptide length': [], 'peptide RMSD': [], 'hamming distance': [], '(pdbid, WT peptide)': []}

for curr__allele_and_length_list, curr__pep_rmsd_list, curr__pdb_indices, conf in [(allele_and_length_list__same_conf, pep_rmsd_list__same_conf, pdb_indices__same_conf, 'same_conf'),
                                                                                   (allele_and_length_list__diff_conf, pep_rmsd_list__diff_conf, pdb_indices__diff_conf, 'diff_conf')]:
                                           
    for (allele, peptide_length), pep_rmsd, pdb_pair in zip(curr__allele_and_length_list, curr__pep_rmsd_list, curr__pdb_indices):

        peptide_length = int(peptide_length)

        pwm_1, pwm_2, pwm_combined, pwm_true, pep_1, pep_2, pdb_1, pdb_2 = get_pwms(final_pwms, df, allele, peptide_length, pdb_pair)

        kld_1 = kl_divergence_of_pwms(pwm_true, pwm_1) / peptide_length
        kld_2 = kl_divergence_of_pwms(pwm_true, pwm_2) / peptide_length
        kld_combined = kl_divergence_of_pwms(pwm_true, pwm_combined) / peptide_length

        kld_absolute_improvement = np.mean([kld_1, kld_2]) - kld_combined

        kld_percent_improvement = ((np.mean([kld_1, kld_2]) - kld_combined) / np.mean([kld_1, kld_2])) * 100

        symmetric_kld_between_the_two_pwms = ((kl_divergence_of_pwms(pwm_1, pwm_2) + kl_divergence_of_pwms(pwm_2, pwm_1)) / 2) / peptide_length

        hdist = hamming_distance(pep_1, pep_2) / peptide_length

        # print((allele, peptide_length, hamming_distance(pep_1, pep_2)))

        all_metrics['kld_absolute_improvement'][conf].append(kld_absolute_improvement)
        all_metrics['kld_percent_improvement'][conf].append(kld_percent_improvement)
        all_metrics['symmetric_kld_between_the_two_pwms'][conf].append(symmetric_kld_between_the_two_pwms)
        all_metrics['hamming_distance_between_the_two_peptides'][conf].append(hdist)
        all_metrics['peptide_length'][conf].append(peptide_length)
        all_metrics['marker'][conf].append(peptide_length_to_marker[peptide_length])
        all_metrics['color'][conf].append(peptide_length_to_color[peptide_length])
        all_metrics['pep_rmsd'][conf].append(pep_rmsd)

        table_for_paper['MHC allele'].append(make_pretty_allele(allele))
        table_for_paper['peptide length'].append(peptide_length)
        table_for_paper['peptide RMSD'].append(pep_rmsd)
        table_for_paper['hamming distance'].append(hamming_distance(pep_1, pep_2))
        table_for_paper['(pdbid, WT peptide)'].append(('; '.join([str((pdb_1.upper(), pep_1)), str((pdb_2.upper(), pep_2))])).replace("'", ""))

table_for_paper = pd.DataFrame(table_for_paper)
table_for_paper.to_csv('table_for_paper.tsv', sep='\t', index=None)

make_distribution_plot(all_metrics['kld_absolute_improvement'], 'KLD/L improvement (bits)', 'plots/pmhc_two_structures_kld_absolute_improvement.png', all_metrics['hamming_distance_between_the_two_peptides'], all_metrics['color'])
make_distribution_plot(all_metrics['kld_percent_improvement'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_kld_percent_improvement.png', all_metrics['hamming_distance_between_the_two_peptides'], all_metrics['color'])
make_distribution_plot(all_metrics['symmetric_kld_between_the_two_pwms'], 'symmetric KLD/L\n between the two PWMs\n(bits)', 'plots/pmhc_two_structures_symmetric_kld_between_the_two_pwms.png', all_metrics['hamming_distance_between_the_two_peptides'], all_metrics['color'])
make_distribution_plot(all_metrics['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', 'plots/pmhc_two_structures_hamming_distance_between_the_two_peptides.png', all_metrics['hamming_distance_between_the_two_peptides'], all_metrics['color'])

make_correlation_plot(all_metrics['kld_absolute_improvement'], 'KLD/L improvement (bits)', 'plots/pmhc_two_structures_kld_absolute_improvement_by_hamming_distance.png', all_metrics['hamming_distance_between_the_two_peptides'], all_metrics['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics['color'])
make_correlation_plot(all_metrics['kld_percent_improvement'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_kld_percent_improvement_by_hamming_distance.png', all_metrics['hamming_distance_between_the_two_peptides'], all_metrics['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics['color'])
make_correlation_plot(all_metrics['symmetric_kld_between_the_two_pwms'], 'symmetric KLD/L\n between the two PWMs\n(bits)', 'plots/pmhc_two_structures_symmetric_kld_between_the_two_pwms_by_hamming_distance.png', all_metrics['hamming_distance_between_the_two_peptides'], all_metrics['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics['color'])

make_correlation_plot(all_metrics['kld_absolute_improvement'], 'KLD/L improvement (bits)', 'plots/pmhc_two_structures_kld_absolute_improvement_by_pep_rmsd.png', all_metrics['pep_rmsd'], all_metrics['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics['color'])
make_correlation_plot(all_metrics['kld_percent_improvement'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_kld_percent_improvement_by_pep_rmsd.png', all_metrics['pep_rmsd'], all_metrics['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics['color'])
make_correlation_plot(all_metrics['symmetric_kld_between_the_two_pwms'], 'symmetric KLD/L\n between the two PWMs\n(bits)', 'plots/pmhc_two_structures_symmetric_kld_between_the_two_pwms_by_pep_rmsd.png', all_metrics['pep_rmsd'], all_metrics['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics['color'])

make_correlation_plot(all_metrics['pep_rmsd'], r'peptide RMSD ($\AA$)', 'plots/pmhc_pep_rmsd_vs_hamming_distance_between_the_two_pwms.png', all_metrics['hamming_distance_between_the_two_peptides'], all_metrics['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics['color'])


# repeat, but low hamming distance
low_hamming_distance_cutoff = 1 / 9

masks = {'same_conf': np.array(all_metrics['hamming_distance_between_the_two_peptides']['same_conf']) == low_hamming_distance_cutoff,
         'diff_conf': np.array(all_metrics['hamming_distance_between_the_two_peptides']['diff_conf']) == low_hamming_distance_cutoff}

all_metrics_low_hamming_distance = {}
for metric in all_metrics:
    all_metrics_low_hamming_distance[metric] = {}
    for conf in all_metrics[metric]:
        all_metrics_low_hamming_distance[metric][conf] = np.array(all_metrics[metric][conf])[masks[conf]]

make_correlation_plot(all_metrics_low_hamming_distance['kld_absolute_improvement'], 'KLD/L improvement (bits)', 'plots/pmhc_two_structures_kld_absolute_improvement_by_hamming_distance__low_hamming_distance.png', all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_low_hamming_distance['color'])
make_correlation_plot(all_metrics_low_hamming_distance['kld_percent_improvement'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_kld_percent_improvement_by_hamming_distance__low_hamming_distance.png', all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_low_hamming_distance['color'])
make_correlation_plot(all_metrics_low_hamming_distance['symmetric_kld_between_the_two_pwms'], 'symmetric KLD/L\n between the two PWMs\n(bits)', 'plots/pmhc_two_structures_symmetric_kld_between_the_two_pwms_by_hamming_distance__low_hamming_distance.png', all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_low_hamming_distance['color'])

make_correlation_plot(all_metrics_low_hamming_distance['kld_absolute_improvement'], 'KLD/L improvement (bits)', 'plots/pmhc_two_structures_kld_absolute_improvement_by_pep_rmsd__low_hamming_distance.png', all_metrics_low_hamming_distance['pep_rmsd'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics_low_hamming_distance['color'])
make_correlation_plot(all_metrics_low_hamming_distance['kld_percent_improvement'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_kld_percent_improvement_by_pep_rmsd__low_hamming_distance.png', all_metrics_low_hamming_distance['pep_rmsd'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics_low_hamming_distance['color'])
make_correlation_plot(all_metrics_low_hamming_distance['symmetric_kld_between_the_two_pwms'], 'symmetric KLD/L\n between the two PWMs\n(bits)', 'plots/pmhc_two_structures_symmetric_kld_between_the_two_pwms_by_pep_rmsd__low_hamming_distance.png', all_metrics_low_hamming_distance['pep_rmsd'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics_low_hamming_distance['color'])
make_correlation_plot(all_metrics_low_hamming_distance['pep_rmsd'], r'peptide RMSD ($\AA$)', 'plots/pmhc_pep_rmsd_vs_hamming_distance_between_the_two_pwms__low_hamming_distance.png', all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_low_hamming_distance['color'])


# repeat, but with high hamming distance
high_hamming_distance_cutoff = 8 / 9

masks = {'same_conf': np.array(all_metrics['hamming_distance_between_the_two_peptides']['same_conf']) == high_hamming_distance_cutoff,
         'diff_conf': np.array(all_metrics['hamming_distance_between_the_two_peptides']['diff_conf']) == high_hamming_distance_cutoff}

all_metrics_high_hamming_distance = {}
for metric in all_metrics:
    all_metrics_high_hamming_distance[metric] = {}
    for conf in all_metrics[metric]:
        all_metrics_high_hamming_distance[metric][conf] = np.array(all_metrics[metric][conf])[masks[conf]]

make_correlation_plot(all_metrics_high_hamming_distance['kld_absolute_improvement'], 'KLD/L improvement (bits)', 'plots/pmhc_two_structures_kld_absolute_improvement_by_hamming_distance__high_hamming_distance.png', all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_high_hamming_distance['color'])
make_correlation_plot(all_metrics_high_hamming_distance['kld_percent_improvement'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_kld_percent_improvement_by_hamming_distance__high_hamming_distance.png', all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_high_hamming_distance['color'])
make_correlation_plot(all_metrics_high_hamming_distance['symmetric_kld_between_the_two_pwms'], 'symmetric KLD/L\n between the two PWMs\n(bits)', 'plots/pmhc_two_structures_symmetric_kld_between_the_two_pwms_by_hamming_distance__high_hamming_distance.png', all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_high_hamming_distance['color'])

make_correlation_plot(all_metrics_high_hamming_distance['kld_absolute_improvement'], 'KLD/L improvement (bits)', 'plots/pmhc_two_structures_kld_absolute_improvement_by_pep_rmsd__high_hamming_distance.png', all_metrics_high_hamming_distance['pep_rmsd'], all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics_high_hamming_distance['color'])
make_correlation_plot(all_metrics_high_hamming_distance['kld_percent_improvement'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_kld_percent_improvement_by_pep_rmsd__high_hamming_distance.png', all_metrics_high_hamming_distance['pep_rmsd'], all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics_high_hamming_distance['color'])
make_correlation_plot(all_metrics_high_hamming_distance['symmetric_kld_between_the_two_pwms'], 'symmetric KLD/L\n between the two PWMs\n(bits)', 'plots/pmhc_two_structures_symmetric_kld_between_the_two_pwms_by_pep_rmsd__high_hamming_distance.png', all_metrics_high_hamming_distance['pep_rmsd'], all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics_high_hamming_distance['color'])
make_correlation_plot(all_metrics_high_hamming_distance['pep_rmsd'], r'peptide RMSD ($\AA$)', 'plots/pmhc_pep_rmsd_vs_hamming_distance_between_the_two_pwms__high_hamming_distance.png', all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_high_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_high_hamming_distance['color'])



# repeat, but low peptide RMSD
low_low_rmsd = 0.05

masks = {'same_conf': np.array(all_metrics['hamming_distance_between_the_two_peptides']['same_conf']) == low_hamming_distance_cutoff,
         'diff_conf': np.array(all_metrics['hamming_distance_between_the_two_peptides']['diff_conf']) == low_hamming_distance_cutoff}

all_metrics_low_hamming_distance = {}
for metric in all_metrics:
    all_metrics_low_hamming_distance[metric] = {}
    for conf in all_metrics[metric]:
        all_metrics_low_hamming_distance[metric][conf] = np.array(all_metrics[metric][conf])[masks[conf]]

make_correlation_plot(all_metrics_low_hamming_distance['kld_absolute_improvement'], 'KLD/L improvement (bits)', 'plots/pmhc_two_structures_kld_absolute_improvement_by_hamming_distance__low_hamming_distance.png', all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_low_hamming_distance['color'])
make_correlation_plot(all_metrics_low_hamming_distance['kld_percent_improvement'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_kld_percent_improvement_by_hamming_distance__low_hamming_distance.png', all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_low_hamming_distance['color'])
make_correlation_plot(all_metrics_low_hamming_distance['symmetric_kld_between_the_two_pwms'], 'symmetric KLD/L\n between the two PWMs\n(bits)', 'plots/pmhc_two_structures_symmetric_kld_between_the_two_pwms_by_hamming_distance__low_hamming_distance.png', all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_low_hamming_distance['color'])

make_correlation_plot(all_metrics_low_hamming_distance['kld_absolute_improvement'], 'KLD/L improvement (bits)', 'plots/pmhc_two_structures_kld_absolute_improvement_by_pep_rmsd__low_hamming_distance.png', all_metrics_low_hamming_distance['pep_rmsd'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics_low_hamming_distance['color'])
make_correlation_plot(all_metrics_low_hamming_distance['kld_percent_improvement'], 'KLD/L improvement (%)', 'plots/pmhc_two_structures_kld_percent_improvement_by_pep_rmsd__low_hamming_distance.png', all_metrics_low_hamming_distance['pep_rmsd'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics_low_hamming_distance['color'])
make_correlation_plot(all_metrics_low_hamming_distance['symmetric_kld_between_the_two_pwms'], 'symmetric KLD/L\n between the two PWMs\n(bits)', 'plots/pmhc_two_structures_symmetric_kld_between_the_two_pwms_by_pep_rmsd__low_hamming_distance.png', all_metrics_low_hamming_distance['pep_rmsd'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], r'peptide RMSD ($\AA$)', all_metrics_low_hamming_distance['color'])
make_correlation_plot(all_metrics_low_hamming_distance['pep_rmsd'], r'peptide RMSD ($\AA$)', 'plots/pmhc_pep_rmsd_vs_hamming_distance_between_the_two_pwms__low_hamming_distance.png', all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], all_metrics_low_hamming_distance['hamming_distance_between_the_two_peptides'], 'hamming distance (%)', all_metrics_low_hamming_distance['color'])


