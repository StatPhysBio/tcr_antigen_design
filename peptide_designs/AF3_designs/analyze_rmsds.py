
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
from scipy.stats import pearsonr, spearmanr

color_list = plt.get_cmap('tab20').colors
blue = color_list[0]
blue_light = color_list[1]
orange = color_list[2]
orange_light = color_list[3]
green = color_list[4]
green_light = color_list[5]
red = color_list[6]
red_light = color_list[7]
purple = color_list[8]
purple_light = color_list[9]
brown = color_list[10]
brown_light = color_list[11]
pink = color_list[12]
pink_light = color_list[13]
gray = color_list[14]
gray_light = color_list[15]
olive = color_list[16]
olive_light = color_list[17]
cyan = color_list[18]
cyan_light = color_list[19]

MODEL_TO_COLOR = {
    ('hermes_py_000', 'fixed'): red,
    ('hermes_py_050', 'fixed'): red_light,
    ('hermes_py_000', 'relaxed'): orange,
    ('hermes_py_050', 'relaxed'): orange_light,
}

MODEL_TO_COLOR_SIMPLER = {
    ('hermes_py_000', 'fixed'): red,
    ('hermes_py_050', 'fixed'): red,
    ('hermes_py_000', 'relaxed'): pink,
    ('hermes_py_050', 'relaxed'): pink,
}

PRETTY_NAME_TO_COLOR = {
    'HERMES-$fixed$': red,
    'HERMES-$relaxed$': pink
}

from matplotlib.patches import Patch

# Create legend handles with LaTeX-style labels and black edges
legend_handles = [
    Patch(facecolor=color, edgecolor='none', label=label)
    for label, color in PRETTY_NAME_TO_COLOR.items()
]

# Dummy plot
fig, ax = plt.subplots()
ax.axis('off')
# Add legend with mathtext rendering
ax.legend(handles=legend_handles, fontsize=14)

plt.savefig('plots/hermes_protocol_legend.png')
plt.savefig('plots/hermes_protocol_legend.pdf')
plt.close()

MODEL_TO_PRETTY_NAME = {
        ('hermes_py_000', 'fixed'): 'HERMES-$fixed$ 0.00',
        ('hermes_py_050', 'fixed'): 'HERMES-$fixed$ 0.50',
        ('hermes_py_000', 'relaxed'): 'HERMES-$relaxed$ 0.00',
        ('hermes_py_050', 'relaxed'): 'HERMES-$relaxed$ 0.50',
}

RMSD_COL_TO_PRETTY_NAME = {
    'pep_rmsd': 'Peptide RMSD',
    'tcr_rmsd': 'TCR RMSD',
}

os.makedirs('plots', exist_ok=True)


df_1 = pd.read_csv('nyeso_mean_rmsds.csv')
df_2 = pd.read_csv('/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/nyeso/post_experiment_analysis/NYES0_test_res__for_paper.csv')
df_2 = df_2.rename(columns={'sequence': 'peptide'})
df_nyeso = pd.merge(df_1, df_2, on='peptide', how='left')


df_1 = pd.read_csv('ebv_mean_rmsds.csv')
df_2 = pd.read_csv('/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/ebv/post_experiment_analysis/EBV_test_res__for_paper.csv')
df_2 = df_2.rename(columns={'sequence': 'peptide'})
df_ebv = pd.merge(df_1, df_2, on='peptide', how='left')


df_1 = pd.read_csv('magea3_titin_mean_rmsds.csv')
df_2 = pd.read_csv('/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/peptide_designs/magea3_and_titin/post_experiment_analysis/MAGE_test_res__for_paper.csv')
df_2 = df_2.rename(columns={'sequence': 'peptide'})
df_magea3 = pd.merge(df_1, df_2, on='peptide', how='left')


def plot_pep_vs_tcr_rmsd(df, name):

    model_and_protocol_list = [('hermes_py_000', 'fixed'),
                                ('hermes_py_050', 'fixed'),
                                ('hermes_py_000', 'relaxed'),
                                ('hermes_py_050', 'relaxed')]

    fontsize = 14

    # Create a scatter plot
    plt.figure(figsize=(2.5, 2.5))
    colors = [MODEL_TO_COLOR_SIMPLER[(model, protocol)] for model, protocol in zip(df['model'].values, df['protocol'].values)]
    for tcr_rmsd, pep_rmsd, color in zip(df['tcr_rmsd'].values, df['pep_rmsd'].values, colors):
        plt.scatter(tcr_rmsd, pep_rmsd, alpha=0.7, facecolor=color, edgecolor='none', s=60)
    plt.xlabel('TCR RMSD', fontsize=fontsize)
    plt.ylabel('Peptide RMSD', fontsize=fontsize)

    # Calculate the correlation coefficient
    r, r_pval = pearsonr(df['tcr_rmsd'], df['pep_rmsd'])
    plt.text(0.9, 0.15, fr'$R$: {r:.2f}', transform=plt.gca().transAxes, fontsize=fontsize, ha='right')
    plt.text(0.95, 0.05, f'p-val: {r_pval:.2f}', transform=plt.gca().transAxes, fontsize=fontsize, ha='right')

    plt.xticks(fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-2)

    plt.grid(ls=':', color='grey', alpha=0.5)

    # Show the plot
    plt.savefig('plots/pep_vs_tcr_rmsd_' + name + '.pdf', bbox_inches='tight')
    plt.savefig('plots/pep_vs_tcr_rmsd_' + name + '.png', bbox_inches='tight')
    plt.close()

def rmsd_boxplot_split_by_model(system_to_df, rmsd_col):

    systems = list(system_to_df.keys())

    model_and_protocol_list = [('hermes_py_000', 'fixed'),
                                ('hermes_py_050', 'fixed'),
                                ('hermes_py_000', 'relaxed'),
                                ('hermes_py_050', 'relaxed')]

    # Plot setup
    fig, ax = plt.subplots(figsize=(6.5, 3))

    # Position tracker
    positions = []
    plot_data = []
    box_colors = []

    # Box width and spacing
    box_width = 0.7
    group_spacing = 0.3  # extra space between groups
    current_position = 1

    fontsize = 14

    # Create mapping from position to raw values for scatter
    scatter_positions = []

    for name in system_to_df:
        df = system_to_df[name]
        for model_and_protocol in model_and_protocol_list:
            data = df.loc[np.logical_and(df['model'] == model_and_protocol[0], df['protocol'] == model_and_protocol[1]), rmsd_col].values
            plot_data.append(data)
            positions.append(current_position)
            box_colors.append(MODEL_TO_COLOR[model_and_protocol])
            scatter_positions.append((current_position, data, box_colors[-1]))
            current_position += 1
        current_position += group_spacing  # Add space after each group
    
    # Create the boxplot
    bp = ax.boxplot(plot_data, positions=positions, widths=box_width, patch_artist=True, showfliers=False)

    # Customize each element
    for patch, color in zip(bp['boxes'], box_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)
        patch.set_linewidth(1)  # Thinner box line

    # Set median line to black
    for median in bp['medians']:
        median.set_color('black')
        median.set_linewidth(1)

    # Thin whiskers and caps
    for line in bp['whiskers'] + bp['caps']:
        line.set_linewidth(1)

    # Overlay scatter points with jitter
    for pos, values, color in scatter_positions:
        jitter = np.random.normal(loc=0, scale=0.1, size=len(values))  # small horizontal jitter
        ax.scatter(pos + jitter, values, alpha=1.0, color=color, edgecolor='none', s=30)

    # Set x-tick labels
    ax.set_xticks([np.mean(positions[i*4:(i+1)*4]) for i in range(len(systems))])
    ax.set_xticklabels(systems, fontsize=fontsize)

    ax.tick_params(axis='y', labelsize=fontsize-2)

    ax.grid(axis='y', ls=':', color='grey', alpha=0.5)

    # Labels and layout
    ax.set_ylabel(RMSD_COL_TO_PRETTY_NAME[rmsd_col], fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(f'plots/rmsd_boxplot_split_by_model_{rmsd_col}.pdf', bbox_inches='tight')
    plt.savefig(f'plots/rmsd_boxplot_split_by_model_{rmsd_col}.png', bbox_inches='tight')
    plt.close()

    # save legend to separate figure
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=color, edgecolor='none', alpha=0.7, label=MODEL_TO_PRETTY_NAME[model]) for model, color in MODEL_TO_COLOR.items()]
    fig_legend, ax_legend = plt.subplots(figsize=(2, 1))
    ax_legend.axis('off')
    ax_legend.legend(handles, [MODEL_TO_PRETTY_NAME[model] for model in MODEL_TO_COLOR], fontsize=fontsize-2)
    plt.savefig(f'plots/rmsd_boxplot_split_by_model_{rmsd_col}_legend.pdf', bbox_inches='tight')
    plt.savefig(f'plots/rmsd_boxplot_split_by_model_{rmsd_col}_legend.png', bbox_inches='tight')
    plt.close(fig_legend)


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


import sys
sys.path.append('/Users/gianmarcovisani/Desktop/tcr_pmhc/tcr_antigen_design/mutation_effects/src')
from global_constants import LOGOMAKER_COLORS

def plot_pwm(pwm, out_path):

    import logomaker
        
    # make the logo show information content, not probability
    information = np.log2(20) + np.sum(pwm * np.log2(pwm + 1e-10), axis=1)

    information_adjusted_pwm = pwm * information[:, np.newaxis]

    information_adjusted_pwm_df = pd.DataFrame(information_adjusted_pwm, index=range(1, pwm.shape[0]+1), columns=[ind_to_ol_size[ind] for ind in range(20)])

    # make the figure
    fontsize = 18
    fig, ax = plt.subplots(figsize=(10, 2))
    ax.set_ylim(0, np.log2(20))
    logomaker.Logo(information_adjusted_pwm_df, ax=ax, color_scheme=LOGOMAKER_COLORS)
    # ax.axis('off')
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.set_ylim([0, np.log2(20)])
    ax.yaxis.set_ticks_position('right')  # Set ticks on the right
    ax.yaxis.set_label_position('right')  # Set the label on the right
    ax.set_ylabel('bits', fontsize=fontsize)
    ax.set_yticks([1, 3])
    ax.tick_params(axis='y', labelsize=fontsize-2)
    ax.set_xticks([])
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()


plot_pep_vs_tcr_rmsd(df_nyeso, 'nyeso')
plot_pep_vs_tcr_rmsd(df_ebv, 'ebv')
plot_pep_vs_tcr_rmsd(df_magea3, 'magea3')


rmsd_boxplot_split_by_model({'NY-ESO': df_nyeso, 'EBV': df_ebv, 'MAGE': df_magea3}, 'pep_rmsd')
rmsd_boxplot_split_by_model({'NY-ESO': df_nyeso, 'EBV': df_ebv, 'MAGE': df_magea3}, 'tcr_rmsd')


# make logoplots of the two conformation groups in NYESO

sequences_nyeso_conf_1 = df_nyeso.loc[df_nyeso['pep_rmsd'] <= 0.8, 'peptide'].values
sequences_nyeso_conf_2 = df_nyeso.loc[df_nyeso['pep_rmsd'] > 0.8, 'peptide'].values

plot_pwm(normalize_pwm(make_pwm_from_sequences(sequences_nyeso_conf_1)), 'plots/pwm_nyeso_conf_pep_rmsd_below_0p8.png')
plot_pwm(normalize_pwm(make_pwm_from_sequences(sequences_nyeso_conf_2)), 'plots/pwm_nyeso_conf_pep_rmsd_above_0p8.png')

print([x.lower() for x in sequences_nyeso_conf_1])
print([x.lower() for x in sequences_nyeso_conf_2])


# make logoplots of the two conformation groups in MAGE

sequences_magea3_conf_1 = df_magea3.loc[df_magea3['tcr_rmsd'] <= 6.0, 'peptide'].values
sequences_magea3_conf_2 = df_magea3.loc[df_magea3['tcr_rmsd'] > 6.0, 'peptide'].values

plot_pwm(normalize_pwm(make_pwm_from_sequences(sequences_magea3_conf_1)), 'plots/pwm_magea3_conf_tcr_rmsd_below_6p0.png')
plot_pwm(normalize_pwm(make_pwm_from_sequences(sequences_magea3_conf_2)), 'plots/pwm_magea3_conf_tcr_rmsd_above_6p0.png')

print([x.lower() for x in sequences_magea3_conf_1])
print([x.lower() for x in sequences_magea3_conf_2])




