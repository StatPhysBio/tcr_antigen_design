
import os
import gzip, pickle
import numpy as np
import pandas as pd
from tqdm import tqdm

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
os.makedirs('plots', exist_ok=True)

# every PWM here follows this numbering
# from hermes.utils.protein_naming import ind_to_ol_size, ol_to_ind_size
ind_to_ol_size = {0: 'G', 1: 'A', 2: 'C', 3: 'S', 4: 'P', 5: 'T', 6: 'V', 7: 'D', 8: 'I', 9: 'L', 10: 'N', 11: 'M', 12: 'Q', 13: 'K', 14: 'E', 15: 'H', 16: 'F', 17: 'R', 18: 'Y', 19: 'W'}
ol_to_ind_size = {ind_to_ol_size[key]: key for key in ind_to_ol_size}

def normalize_pwm(pwm):
    return pwm / np.sum(pwm, axis=1, keepdims=True)

def make_pwm_from_sequences(sequences):
    L = len(sequences[0])
    pwm = np.zeros((L, 20))

    for sequence in sequences:
        assert len(sequence) == L, 'All sequences must have the same length'
        if 'X' in sequence:
            print(f'Warning: found X in sequence {sequence}. Skipping.')
            continue
        for seq_i in range(L):
            aa = sequence[seq_i]
            aa_i = ol_to_ind_size[aa]
            pwm[seq_i, aa_i] += 1
    
    return normalize_pwm(pwm)


# with gzip.open('/gscratch/spe/gvisan01/peptide_mhc/mhc_motif_atlas/mhc_motif_pwms.pkl.gz', 'rb') as f:
#     mhc_motif_pwms = pickle.load(f)

# df = pd.read_csv('pmhc_class_1_crystal_structures.csv')

# def fix_allele(allele):
    
#     if isinstance(allele, str):
#         if '*' in allele:
#             return allele[:1] + ''.join(allele[2:].split(':')[:2])
#         else:
#             return allele[:2] + '-' + allele[2:]
#     else:
#         return allele

# df['mhc_allele_final'] = df['mhc_allele_longer'].apply(fix_allele)
# df['peptide_length'] = df['wt_peptide'].apply(len)

# df = df.dropna(subset=['mhc_allele_final'])

# # Exclude rows where wt_peptide contains an X. I don't know what it means.
# pep_does_not_contains_x = df['wt_peptide'].apply(lambda pep: 'X' not in pep)
# df = df[pep_does_not_contains_x]


# print(len(df))

# final_pwms = {}

# groups = df.groupby(['mhc_allele_final', 'peptide_length'])

# for group, group_df in tqdm(groups, total=len(groups)):
#     allele, pep_length = group

#     if allele not in final_pwms:
#         final_pwms[allele] = {}
#     if pep_length not in final_pwms[allele]:
#         final_pwms[allele][pep_length] = {}
    
#     try:
#         mhc_motif_pwm = normalize_pwm(mhc_motif_pwms['class_I'][allele][pep_length])
#     except KeyError:
#         print(f'Warning: Could not find MHC motif PWM for {group}. Skipping group of size {len(group_df)}.')
#         continue

#     # now making hermes' pwm
#     pwms = []
#     rows_to_keep = []
#     for i_row, row in group_df.iterrows():
#         pdb = row['pdbid']
#         wt_pep = row['wt_peptide']
#         old_allele = row['mhc_allele']
#         try:
#             pwm = normalize_pwm(pd.read_csv(f'pwm_csv_files/mhc_crystal_hcnn_fixed_structure/hermes_py_000/{pdb}__{wt_pep}__{old_allele}.csv', index_col=0).values)
#         except FileNotFoundError:
#             print(f'Warning: {pdb}__{wt_pep}__{old_allele}.csv not found. Hermes likely failed on it.')
#             continue

#         if pwm.shape[0] != pep_length:
#             print('Warning: shape mismatch in PWM, something likely wrong with structure or structure parsing.')
#             continue
#         pwms.append(pwm)
#         rows_to_keep.append(i_row)

#     filtered_group_df = group_df.loc[np.array(rows_to_keep)]
    
#     wt_pep_in_struc_pwm = make_pwm_from_sequences(filtered_group_df['wt_peptide'].values)

#     try:
#         hermes_pwm = np.mean(np.stack(pwms, axis=0), axis=0) # it's already normalized
#     except:
#         print(group)
#         print(len(pwms))
#         for pwm in pwms:
#             print('\t', pwm.shape)
#         raise

#     final_pwms[allele][pep_length] = {
#         'mhc_motif_pwm': mhc_motif_pwm,
#         'wt_pep_in_struc_pwm': wt_pep_in_struc_pwm,
#         'hermes_py_000_pwm': hermes_pwm,
#         'num_structures': len(filtered_group_df)
#     }

# with gzip.open('pmhc_class_1_crystal_structures_pwms.pkl.gz', 'wb') as f:
#     pickle.dump(final_pwms, f)


## make three barplots, in one figure, on top of each other
## on the x-axis are the alleles. ordered by human and mouse, and then by number of structures within (or alphabetically?)
## on the three y-axes show:
##      1) number of structures
##      2) entropy of the three pwms (or entropy difference from the MHC Motif PWM? that's probably best)
##      3) relative entropy (kl-divergence) to MHC Motif PWM

## also, make logoplots for each allele, showing information content


with gzip.open('pmhc_class_1_crystal_structures_pwms.pkl.gz', 'rb') as f:
    final_pwms = pickle.load(f)

outdir = 'pmhc_pwms'
os.makedirs(outdir, exist_ok=True)

def entropy_of_pwm(pwm, epsilon=1e-12):
    # Clip to avoid log(0)
    prob_clipped = np.clip(pwm, epsilon, 1.0)
    # entropies = -np.sum(prob_clipped * np.log2(prob_clipped), axis=1)  # shape: (N,)
    # return np.mean(entropies)
    return -np.sum(prob_clipped * np.log2(prob_clipped))

def kl_divergence_of_pwms(P, Q, epsilon=1e-12):
    P_clipped = np.clip(P, epsilon, 1.0)
    Q_clipped = np.clip(Q, epsilon, 1.0)
    
    # kl_div = np.sum(P_clipped * np.log2(P_clipped / Q_clipped), axis=1)  # shape: (N,)
    # return np.mean(kl_div)
    return np.sum(P_clipped * np.log2(P_clipped / Q_clipped))

def plot_pwm(pwm, out_path):

    os.makedirs(outdir, exist_ok=True)

    import logomaker
        
    # make the logo show information content, not probability
    information = np.log2(20) + np.sum(pwm * np.log2(pwm + 1e-10), axis=1)

    information_adjusted_pwm = pwm * information[:, np.newaxis]

    information_adjusted_pwm_df = pd.DataFrame(information_adjusted_pwm, index=range(1, pwm.shape[0]+1), columns=[ind_to_ol_size[ind] for ind in range(20)])

    
    # make the figure
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.set_ylim(0, np.log2(20))
    logomaker.Logo(information_adjusted_pwm_df, ax=ax)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


# pre-compute values entropies
allele_and_length_list = []
num_structures_list = []
mhc_motif_entropy_list = []
wt_pep_in_struc_entropy_list = []
hermes_py_000_entropy_list = []
hermes_py_000_to_mhc_motif_rel_entropy_list = []
wt_pep_in_struc_to_mhc_motif_rel_entropy_list = []

for allele in tqdm(final_pwms):
    for pep_length in final_pwms[allele]:
        
        if len(final_pwms[allele][pep_length]) == 0: # empty dict
            continue

        allele_and_length_list.append((allele, pep_length))

        num_structures = final_pwms[allele][pep_length]['num_structures']

        num_structures_list.append(num_structures)
        mhc_motif_pwm = final_pwms[allele][pep_length]['mhc_motif_pwm']
        wt_pep_in_struc_pwm = final_pwms[allele][pep_length]['wt_pep_in_struc_pwm']
        hermes_py_000_pwm = final_pwms[allele][pep_length]['hermes_py_000_pwm']

        mhc_motif_entropy_list.append(entropy_of_pwm(mhc_motif_pwm))
        wt_pep_in_struc_entropy_list.append(entropy_of_pwm(wt_pep_in_struc_pwm))
        hermes_py_000_entropy_list.append(entropy_of_pwm(hermes_py_000_pwm))

        wt_pep_in_struc_to_mhc_motif_rel_entropy_list.append(kl_divergence_of_pwms(wt_pep_in_struc_pwm, mhc_motif_pwm))
        hermes_py_000_to_mhc_motif_rel_entropy_list.append(kl_divergence_of_pwms(hermes_py_000_pwm, mhc_motif_pwm))

        # plot pwms!
        os.makedirs(f'pmhc_pwms/mhc_motif', exist_ok=True)
        os.makedirs(f'pmhc_pwms/wt_pep_in_struc', exist_ok=True)
        os.makedirs(f'pmhc_pwms/hermes_py_000', exist_ok=True)
        
        plot_pwm(mhc_motif_pwm, f'pmhc_pwms/mhc_motif/mhc_motif__{allele}__{pep_length}.png')
        plot_pwm(wt_pep_in_struc_pwm, f'pmhc_pwms/wt_pep_in_struc/wt_pep_in_struc__{allele}__{pep_length}.png')
        plot_pwm(hermes_py_000_pwm, f'pmhc_pwms/hermes_py_000/hermes_py_000__{allele}__{pep_length}.png')


allele_and_length_list = np.array(allele_and_length_list)
num_structures_list = np.array(num_structures_list)
mhc_motif_entropy_list = np.array(mhc_motif_entropy_list)
wt_pep_in_struc_entropy_list = np.array(wt_pep_in_struc_entropy_list)
hermes_py_000_entropy_list = np.array(hermes_py_000_entropy_list)
wt_pep_in_struc_to_mhc_motif_rel_entropy_list = np.array(wt_pep_in_struc_to_mhc_motif_rel_entropy_list)
hermes_py_000_to_mhc_motif_rel_entropy_list = np.array(hermes_py_000_to_mhc_motif_rel_entropy_list)

## have at least k structures
mask = num_structures_list >= 2

allele_and_length_list = allele_and_length_list[mask]
num_structures_list = num_structures_list[mask]
mhc_motif_entropy_list = mhc_motif_entropy_list[mask]
wt_pep_in_struc_entropy_list = wt_pep_in_struc_entropy_list[mask]
hermes_py_000_entropy_list = hermes_py_000_entropy_list[mask]
wt_pep_in_struc_to_mhc_motif_rel_entropy_list = wt_pep_in_struc_to_mhc_motif_rel_entropy_list[mask]
hermes_py_000_to_mhc_motif_rel_entropy_list = hermes_py_000_to_mhc_motif_rel_entropy_list[mask]

print(mhc_motif_entropy_list)
print(wt_pep_in_struc_entropy_list)
print(hermes_py_000_entropy_list)



def make_pretty_allele_and_length_name(allele_and_length):
    allele, length = allele_and_length
    if '-' in allele: # mouse
        pretty_allele = allele
    else:
        pretty_allele = allele[0] + '*' + allele[1:3] + ':' + allele[3:]
    return f'{pretty_allele}, L={length}'



fontsize = 18

global_bar_width = 0.75

ncols = 1
nrows = 4
colsize = len(allele_and_length_list) * 0.5
rowsize = 3.5
fig, axs = plt.subplots(figsize=(ncols*colsize, nrows*rowsize), ncols=ncols, nrows=nrows, sharex=True, sharey=False)

x = np.arange(len(allele_and_length_list))

# num structures
ax = axs[0]
ax.bar(x, num_structures_list, width=global_bar_width, color='grey')
ax.grid(axis='both', ls='--', alpha=0.5)
ax.set_ylabel('Number of\npMHC structures', fontsize=fontsize)
ax.tick_params(axis='y', labelsize=fontsize-2)

# entropy difference from mhc motif pwm
ax = axs[1]

num_groups = len(allele_and_length_list)
num_bars = 2
bar_width = global_bar_width / num_bars  # Automatically spread bars within group
x = np.arange(num_groups)  # Base x-locations

ax.set_ylabel('Relative entropy\nto MHC Motif PWM', fontsize=fontsize)
entropy_differences = [
    wt_pep_in_struc_to_mhc_motif_rel_entropy_list,
    hermes_py_000_to_mhc_motif_rel_entropy_list
]

colors = ['blue', 'orange']

for i in range(num_bars):
    offset = (i - num_bars / 2) * bar_width + bar_width / 2
    ax.bar(x + offset, entropy_differences[i], width=bar_width, color=colors[i])

# ax.axhline(0.0, ls='--', color='black')
ax.grid(axis='both', ls='--', alpha=0.5)
ax.tick_params(axis='y', labelsize=fontsize-2)


# entropy difference from mhc motif pwm
ax = axs[2]

num_groups = len(allele_and_length_list)
num_bars = 2
bar_width = global_bar_width / num_bars  # Automatically spread bars within group
x = np.arange(num_groups)  # Base x-locations

ax.set_ylabel('Entropy difference\nfrom MHC Motif PWM', fontsize=fontsize)
entropy_differences = [
    wt_pep_in_struc_entropy_list - mhc_motif_entropy_list,
    hermes_py_000_entropy_list - mhc_motif_entropy_list
]

colors = ['blue', 'orange']
labels = ['PWM from WT peptides', 'PWM from HERMES-fixed']

for i in range(num_bars):
    offset = (i - num_bars / 2) * bar_width + bar_width / 2
    ax.bar(x + offset, entropy_differences[i], width=bar_width, color=colors[i], label=labels[i])

# ax.axhline(0.0, ls='--', color='black')
ax.grid(axis='both', ls='--', alpha=0.5)


ax.set_xticks(x)
ax.set_xticklabels(map(make_pretty_allele_and_length_name, allele_and_length_list), rotation=70, ha='right')
ax.tick_params(axis='y', labelsize=fontsize-2)
ax.tick_params(axis='x', labelsize=fontsize-2)

# ax.legend(loc='lower right', fontsize=fontsize-1)

# just entropies!
ax = axs[3]

num_groups = len(allele_and_length_list)
num_bars = 3
bar_width = global_bar_width / num_bars  # Automatically spread bars within group
x = np.arange(num_groups)  # Base x-locations

ax.set_ylabel('Entropy', fontsize=fontsize)
entropies = [
    wt_pep_in_struc_entropy_list,
    hermes_py_000_entropy_list,
    mhc_motif_entropy_list
]

colors = ['blue', 'orange', 'green']
labels = ['PWM from WT peptides', 'PWM from HERMES-fixed', 'PWM from MHC Motif Atlas']

for i in range(num_bars):
    offset = (i - num_bars / 2) * bar_width + bar_width / 2
    ax.bar(x + offset, entropies[i], width=bar_width, color=colors[i], label=labels[i])

# ax.axhline(0.0, ls='--', color='black')
ax.grid(axis='both', ls='--', alpha=0.5)


ax.set_xticks(x)
ax.set_xticklabels(map(make_pretty_allele_and_length_name, allele_and_length_list), rotation=70, ha='right')
ax.tick_params(axis='y', labelsize=fontsize-1)
ax.tick_params(axis='x', labelsize=fontsize-1)

ax.legend(loc='upper right', fontsize=fontsize-1)

plt.tight_layout()
plt.savefig('pmhc_pwm_entropy_comparisons.png')
plt.savefig('pmhc_pwm_entropy_comparisons.pdf')
plt.close()




