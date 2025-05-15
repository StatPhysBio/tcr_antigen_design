
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def plot_distribution(scores, wt_score, outpath, peptide_length):

    fig, ax = plt.subplots(figsize=(5, 4))
    ax.hist(scores, bins=50, color='blue', alpha=0.7)
    ax.axvline(wt_score, color='red', linewidth=2)
    ax.set_xlabel(r'$E_{pep}/$' + f'${str(peptide_length)}$')
    ax.set_ylabel('Count')
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close(fig)


## nyeso

df = pd.read_csv('../nyeso/nyeso_full_copy/')

