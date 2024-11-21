
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

if __name__ == '__main__':

    tsvfiles = ['wildtype/wildtype_w_pae.tsv',
                # 'blosum62/sample_peptides_from_blosum62__temperature=1.0_w_pae.tsv',
                # 'blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae.tsv',
                'hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae.tsv',
                'hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae.tsv',
                'hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae.tsv',
                'hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae.tsv',
                'mhc_pwm/mhc_motif_peptides_w_pae.tsv',
                'proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7_w_pae.tsv',
                'proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7_w_pae.tsv']

    neg_pae_scores_original_design = []
    neg_pae_scores_designs_with_e1 = []

    for tsvfile in tsvfiles:
        df = pd.read_csv(tsvfile, sep='\t')
        neg_pae_scores_designs_with_e1.extend(list(-df['pmhc_tcr_pae'].values))

        df = pd.read_csv(os.path.join('../magea3/', tsvfile), sep='\t')
        neg_pae_scores_original_design.extend(list(-df['pmhc_tcr_pae'].values))
    
    ## make square scatter plot, with xlim and ylim equal, diagonal line, neg_pae_scores_designs_with_e1 vs. neg_pae_scores_original_design
    plt.figure(figsize=(5, 5))
    plt.scatter(neg_pae_scores_original_design, neg_pae_scores_designs_with_e1, alpha=0.5, s=10)

    ax = plt.gca()
    # set xlim equal to ylim, and plot diagonal line
    xlim = ax.get_xlim()[0], ax.get_xlim()[1]
    ylim = ax.get_ylim()[0], ax.get_ylim()[1]
    ax.set_xlim(min(xlim[0], ylim[0]), max(xlim[1], ylim[1]))
    ax.set_ylim(min(xlim[0], ylim[0]), max(xlim[1], ylim[1]))
    ax.plot([0,1],[0,1], c='k', transform=ax.transAxes)

    plt.tick_params(axis='both', which='major', labelsize=14)
    plt.xlabel(f'Original Designs\n(mean = {np.nanmean(neg_pae_scores_original_design):.2f}; std = {np.nanstd(neg_pae_scores_original_design):.2f})', fontsize=16)
    plt.ylabel(f'Designs with E1\n(mean = {np.nanmean(neg_pae_scores_designs_with_e1):.2f}; std = {np.nanstd(neg_pae_scores_designs_with_e1):.2f})', fontsize=16)
    plt.title('Negative TCRDock PAE Scores', fontsize=16)

    plt.tight_layout()
    plt.savefig('neg_pae_scores__designs_with_e1_vs__original_designs.png')
    plt.close()


