
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42

import logomaker

sys.path.append('../../../mutation_effects/src')
from global_constants import LOGOMAKER_COLORS

# every PWM here follows this numbering
# from hermes.utils.protein_naming import ind_to_ol_size, ol_to_ind_size
ind_to_ol_size = {0: 'G', 1: 'A', 2: 'C', 3: 'S', 4: 'P', 5: 'T', 6: 'V', 7: 'D', 8: 'I', 9: 'L', 10: 'N', 11: 'M', 12: 'Q', 13: 'K', 14: 'E', 15: 'H', 16: 'F', 17: 'R', 18: 'Y', 19: 'W'}
ol_to_ind_size = {ind_to_ol_size[key]: key for key in ind_to_ol_size}

def normalize_pwm(pwm):
    return pwm / np.sum(pwm, axis=1, keepdims=True)

def plot_pwm(pwm, out_path):
        
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


workingdir = 'hermes_py_000'

for file in os.listdir(workingdir):
    if file.endswith('.csv'):
        pwm = normalize_pwm(pd.read_csv(os.path.join(workingdir, file), index_col=0).values)
        plot_pwm(pwm, os.path.join(workingdir, file.replace('.csv', '.png')))

