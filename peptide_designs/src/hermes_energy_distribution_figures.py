
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


def plot_distribution(scores, wt_scores, outpath, peptide_length):
    fontsize = 20
    fig, ax = plt.subplots(figsize=(4, 3))
    ax.hist(scores / peptide_length, bins=50, color='blue', alpha=0.7)
    for wt_score in wt_scores:
        ax.axvline(wt_score / peptide_length, color='red', linewidth=2)
    ax.set_xlabel(r'$E_{pep}/$' + f'${str(peptide_length)}$', fontsize=fontsize)
    ax.set_ylabel('number of\n designs', fontsize=fontsize)
    ax.tick_params(axis='both', labelsize=fontsize-1)
    plt.tight_layout()
    plt.savefig(outpath)
    plt.close()


## HERMES-fixed

## nyeso

df_hermes_fixed_000 = pd.read_csv('../nyeso/nyeso_full_copy/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t')
df_hermes_fixed_050 = pd.read_csv('../nyeso/nyeso_full_copy/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')
df_hermes_fixed_000_wt = pd.read_csv('../nyeso/nyeso_full_copy/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t')
df_hermes_fixed_050_wt = pd.read_csv('../nyeso/nyeso_full_copy/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_050.tsv', sep='\t')

plot_distribution(df_hermes_fixed_000['pnE'].values,
                  df_hermes_fixed_000_wt['pnE'].values,
                  '../nyeso/hermes_py_000_fixed_scores_distribution.png',
                  len(df_hermes_fixed_000_wt['peptide'].values[0]))

plot_distribution(df_hermes_fixed_050['pnE'].values,
                  df_hermes_fixed_050_wt['pnE'].values,
                  '../nyeso/hermes_py_050_fixed_scores_distribution.png',
                  len(df_hermes_fixed_050_wt['peptide'].values[0]))


## ebv

df_hermes_fixed_000 = pd.concat([pd.read_csv('../ebv/hpvg/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                 pd.read_csv('../ebv/hpvg_q5/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t')])
df_hermes_fixed_050 = pd.concat([pd.read_csv('../ebv/hpvg/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t'),
                                 pd.read_csv('../ebv/hpvg_q5/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])
df_hermes_fixed_000_wt = pd.concat([pd.read_csv('../ebv/hpvg/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t'),
                                    pd.read_csv('../ebv/hpvg_q5/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t')])
df_hermes_fixed_050_wt = pd.concat([pd.read_csv('../ebv/hpvg/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_050.csv', sep=','),
                                    pd.read_csv('../ebv/hpvg_q5/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_050.tsv', sep='\t')])

plot_distribution(df_hermes_fixed_000['pnE'].values,
                  df_hermes_fixed_000_wt['pnE'].values,
                  '../ebv/hermes_py_000_fixed_scores_distribution.png',
                  len(df_hermes_fixed_000_wt['peptide'].values[0]))

plot_distribution(df_hermes_fixed_050['pnE'].values,
                  df_hermes_fixed_050_wt['pnE'].values,
                  '../ebv/hermes_py_050_fixed_scores_distribution.png',
                  len(df_hermes_fixed_050_wt['peptide'].values[0]))


## magea3

df_hermes_fixed_000 = pd.concat([pd.read_csv('../magea3_and_titin/magea3/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                 pd.read_csv('../magea3_and_titin/magea3_with_e1/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                 pd.read_csv('../magea3_and_titin/titin/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t')])

df_hermes_fixed_050 = pd.concat([pd.read_csv('../magea3_and_titin/magea3/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t'),
                                 pd.read_csv('../magea3_and_titin/magea3_with_e1/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t'),
                                 pd.read_csv('../magea3_and_titin/titin/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])

df_hermes_fixed_000_wt = pd.concat([pd.read_csv('../magea3_and_titin/magea3/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t'),
                                    pd.read_csv('../magea3_and_titin/magea3_with_e1/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t'),
                                    pd.read_csv('../magea3_and_titin/titin/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t')])

df_hermes_fixed_050_wt = pd.concat([pd.read_csv('../magea3_and_titin/magea3/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_050.tsv', sep='\t'),
                                    pd.read_csv('../magea3_and_titin/magea3_with_e1/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_050.tsv', sep='\t'),
                                    pd.read_csv('../magea3_and_titin/titin/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_050.tsv', sep='\t')])

plot_distribution(df_hermes_fixed_000['pnE'].values,
                  df_hermes_fixed_000_wt['pnE'].values,
                  '../magea3_and_titin/hermes_py_000_fixed_scores_distribution.png',
                  len(df_hermes_fixed_000_wt['peptide'].values[0]))

plot_distribution(df_hermes_fixed_050['pnE'].values,
                  df_hermes_fixed_050_wt['pnE'].values,
                  '../magea3_and_titin/hermes_py_050_fixed_scores_distribution.png',
                  len(df_hermes_fixed_050_wt['peptide'].values[0]))



## HERMES-relaxed

def get_wt_scores(df_with_scores, wt_sequences, wt_pdbs):
    wt_scores = []
    for seq, pdb in zip(wt_sequences, wt_pdbs):
        wt_scores.append(df_with_scores.loc[np.logical_and(df_with_scores['sequence'] == seq, df_with_scores['pdb'] == pdb)]['pnE'].values[0])
    return wt_scores

## nyeso

df_hermes_relaxed_000 = pd.read_csv('../nyeso/nyeso_full_copy/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t')
df_hermes_relaxed_050 = pd.read_csv('../nyeso/nyeso_full_copy/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')

# take the wt scores from where we already computed them
wt_sequences = pd.read_csv('../nyeso/nyeso_full_copy/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t')['peptide'].values
wt_pdbs = pd.read_csv('../nyeso/nyeso_full_copy/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t')['pdbid'].values
df_wt_000_temp = pd.read_csv('../nyeso/hermes_scores/hermes_py_000/selected_designs_for_hermes_scoring__hermes_py_000__relaxed.csv')
df_wt_050_temp = pd.read_csv('../nyeso/hermes_scores/hermes_py_050/selected_designs_for_hermes_scoring__hermes_py_050__relaxed.csv')
hermes_relaxed_000_wt_scores = get_wt_scores(df_wt_000_temp, wt_sequences, wt_pdbs)
hermes_relaxed_050_wt_scores = get_wt_scores(df_wt_050_temp, wt_sequences, wt_pdbs)


plot_distribution(df_hermes_relaxed_000['pnE'].values,
                  hermes_relaxed_000_wt_scores,
                  '../nyeso/hermes_py_000_relaxed_scores_distribution.png',
                  len(wt_sequences[0]))

plot_distribution(df_hermes_relaxed_050['pnE'].values,
                  hermes_relaxed_050_wt_scores,
                  '../nyeso/hermes_py_050_relaxed_scores_distribution.png',
                  len(wt_sequences[0]))


## ebv

df_hermes_relaxed_000 = pd.concat([pd.read_csv('../ebv/hpvg/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                   pd.read_csv('../ebv/hpvg_q5/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t')])
df_hermes_relaxed_050 = pd.concat([pd.read_csv('../ebv/hpvg/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t'),
                                   pd.read_csv('../ebv/hpvg_q5/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])

# take the wt scores from where we already computed them
wt_sequences = pd.concat([pd.read_csv('../ebv/hpvg/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t'),
                          pd.read_csv('../ebv/hpvg_q5/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t')])['peptide'].values
wt_pdbs = pd.concat([pd.read_csv('../ebv/hpvg/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t'),
                          pd.read_csv('../ebv/hpvg_q5/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t')])['pdbid'].values
df_wt_000_temp = pd.read_csv('../ebv/hermes_scores/hermes_py_000/selected_designs_for_hermes_scoring__hermes_py_000__relaxed.csv')
df_wt_050_temp = pd.read_csv('../ebv/hermes_scores/hermes_py_050/selected_designs_for_hermes_scoring__hermes_py_050__relaxed.csv')
hermes_relaxed_000_wt_scores = get_wt_scores(df_wt_000_temp, wt_sequences, wt_pdbs)
hermes_relaxed_050_wt_scores = get_wt_scores(df_wt_050_temp, wt_sequences, wt_pdbs)


plot_distribution(df_hermes_relaxed_000['pnE'].values,
                  hermes_relaxed_000_wt_scores,
                  '../ebv/hermes_py_000_relaxed_scores_distribution.png',
                  len(wt_sequences[0]))

plot_distribution(df_hermes_relaxed_050['pnE'].values,
                  hermes_relaxed_050_wt_scores,
                  '../ebv/hermes_py_050_relaxed_scores_distribution.png',
                  len(wt_sequences[0]))


## magea3

df_hermes_relaxed_000 = pd.concat([pd.read_csv('../magea3_and_titin/magea3/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                   pd.read_csv('../magea3_and_titin/magea3_with_e1/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                   pd.read_csv('../magea3_and_titin/titin/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t')])
df_hermes_relaxed_050 = pd.concat([pd.read_csv('../magea3_and_titin/magea3/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t'),
                                   pd.read_csv('../magea3_and_titin/magea3_with_e1/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t'),
                                   pd.read_csv('../magea3_and_titin/titin/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])

# take the wt scores from where we already computed them
wt_sequences = pd.concat([pd.read_csv('../magea3_and_titin/magea3/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t'),
                          pd.read_csv('../magea3_and_titin/magea3_with_e1/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t'),
                          pd.read_csv('../magea3_and_titin/titin/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t')])['peptide'].values
wt_pdbs = pd.concat([pd.read_csv('../magea3_and_titin/magea3/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t'),
                          pd.read_csv('../magea3_and_titin/magea3_with_e1/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t'),
                          pd.read_csv('../magea3_and_titin/titin/wildtype/wildtype_w_pae_w_blosum_w_hermes_fixed_000.tsv', sep='\t')])['pdbid'].values
df_wt_000_temp = pd.read_csv('../magea3_and_titin/hermes_scores/hermes_py_000/selected_designs_for_hermes_scoring__hermes_py_000__relaxed.csv')
df_wt_050_temp = pd.read_csv('../magea3_and_titin/hermes_scores/hermes_py_050/selected_designs_for_hermes_scoring__hermes_py_050__relaxed.csv')
hermes_relaxed_000_wt_scores = get_wt_scores(df_wt_000_temp, wt_sequences, wt_pdbs)
hermes_relaxed_050_wt_scores = get_wt_scores(df_wt_050_temp, wt_sequences, wt_pdbs)


plot_distribution(df_hermes_relaxed_000['pnE'].values,
                  hermes_relaxed_000_wt_scores,
                  '../magea3_and_titin/hermes_py_000_relaxed_scores_distribution.png',
                  len(wt_sequences[0]))

plot_distribution(df_hermes_relaxed_050['pnE'].values,
                  hermes_relaxed_050_wt_scores,
                  '../magea3_and_titin/hermes_py_050_relaxed_scores_distribution.png',
                  len(wt_sequences[0]))
