

'''

96 total
1 MageA3
1 Titin

Below: sample from peptides below PAE_CUTOFF TCRdock PAE, and is binder by NetMHCPan, mix No-noise and Noise designs

16 from fixed structure MageA3
16 from simulated annealing Mage A3
16 from fixed structure Titin
16 from simulated annealing Titin
15 from fixed structure MageA3 with forced E1
15 from simulated annealing MAgeA3 with forced E1

Make plot of Hamming distance from both Mage and Titin of all the designs

Make plot of negative tcrdock pae of all the designs

'''

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as cm

def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    return sum([seq1[i] != seq2[i] for i in range(len(seq1))])

MAGE_SEQ = 'EVDPIGHLY'
TITIN_SEQ = 'ESDPIVAQY'

SEED = 43
np.random.seed(SEED)

PAE_CUTOFF = 5.2


if __name__ == '__main__':

    df_mage_wt = pd.read_csv('magea3/wildtype/wildtype_w_pae_w_blosum.tsv', sep='\t')
    df_titin_wt = pd.read_csv('titin/wildtype/wildtype_w_pae_w_blosum.tsv', sep='\t')

    df_mage_fixed_structure = pd.concat([pd.read_csv('magea3/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                         pd.read_csv('magea3/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])

    df_mage_sim_anneal = pd.concat([pd.read_csv('magea3/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                         pd.read_csv('magea3/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])

    df_titin_fixed_structure = pd.concat([pd.read_csv('titin/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                     pd.read_csv('titin/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])

    df_titin_sim_anneal = pd.concat([pd.read_csv('titin/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                     pd.read_csv('titin/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])

    df_mage_with_e1_fixed_structure = pd.concat([pd.read_csv('magea3_with_e1/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                                 pd.read_csv('magea3_with_e1/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])

    df_mage_with_e1_sim_anneal = pd.concat([pd.read_csv('magea3_with_e1/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                            pd.read_csv('magea3_with_e1/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])

    all_sequences = []
    all_paes = []
    all_model_types = []

    sequences_wildtype = [MAGE_SEQ, TITIN_SEQ]
    paes_wildtype = [df_mage_wt['pmhc_tcr_pae'][0], df_titin_wt['pmhc_tcr_pae'][0]]
    all_sequences.extend(sequences_wildtype)
    all_paes.extend(paes_wildtype)
    all_model_types.extend(['magea3_wt', 'titin_wt'])
    
    df_mage_fixed_structure = df_mage_fixed_structure.loc[np.logical_and(df_mage_fixed_structure['is_binder_by_netmhc_pan'] == 1, df_mage_fixed_structure['pmhc_tcr_pae'] < PAE_CUTOFF)]
    min_hamming_distances = np.array([min([hamming_distance(seq, WT_SEQ) for WT_SEQ in sequences_wildtype]) for seq in df_mage_fixed_structure['peptide'].values])
    df_mage_fixed_structure = df_mage_fixed_structure.loc[min_hamming_distances >= 3]
    indices = np.random.choice(len(df_mage_fixed_structure), size=16, replace=False)
    sequences_mage_fixed_structure = list(df_mage_fixed_structure['peptide'].values[indices])
    while len(set(sequences_mage_fixed_structure).intersection(set(all_sequences))) > 0 or len(set(sequences_mage_fixed_structure)) < len(sequences_mage_fixed_structure):
        indices = np.random.choice(len(df_mage_fixed_structure), size=16, replace=False)
        sequences_mage_fixed_structure = list(df_mage_fixed_structure['peptide'].values[indices])
    paes_mage_fixed_structure = list(df_mage_fixed_structure['pmhc_tcr_pae'].values[indices])
    all_sequences.extend(sequences_mage_fixed_structure)
    all_paes.extend(paes_mage_fixed_structure)
    all_model_types.extend(list(df_mage_fixed_structure['hcnn_model'].values[indices]))

    df_mage_sim_anneal = df_mage_sim_anneal.loc[np.logical_and(df_mage_sim_anneal['is_binder_by_netmhc_pan'] == 1, df_mage_sim_anneal['pmhc_tcr_pae'] < PAE_CUTOFF)]
    min_hamming_distances = np.array([min([hamming_distance(seq, WT_SEQ) for WT_SEQ in sequences_wildtype]) for seq in df_mage_sim_anneal['peptide'].values])
    df_mage_sim_anneal = df_mage_sim_anneal.loc[min_hamming_distances >= 3]
    indices = np.random.choice(len(df_mage_sim_anneal), size=16, replace=False)
    sequences_mage_sim_anneal = list(df_mage_sim_anneal['peptide'].values[indices])
    while len(set(sequences_mage_sim_anneal).intersection(set(all_sequences))) > 0 or len(set(sequences_mage_sim_anneal)) < len(sequences_mage_sim_anneal):
        indices = np.random.choice(len(df_mage_sim_anneal), size=16, replace=False)
        sequences_mage_sim_anneal = list(df_mage_sim_anneal['peptide'].values[indices])
    paes_mage_sim_anneal = list(df_mage_sim_anneal['pmhc_tcr_pae'].values[indices])
    all_sequences.extend(sequences_mage_sim_anneal)
    all_paes.extend(paes_mage_sim_anneal)
    all_model_types.extend(list(df_mage_sim_anneal['hcnn_model'].values[indices]))

    df_titin_fixed_structure = df_titin_fixed_structure.loc[np.logical_and(df_titin_fixed_structure['is_binder_by_netmhc_pan'] == 1, df_titin_fixed_structure['pmhc_tcr_pae'] < PAE_CUTOFF)]
    min_hamming_distances = np.array([min([hamming_distance(seq, WT_SEQ) for WT_SEQ in sequences_wildtype]) for seq in df_titin_fixed_structure['peptide'].values])
    df_titin_fixed_structure = df_titin_fixed_structure.loc[min_hamming_distances >= 3]
    indices = np.random.choice(len(df_titin_fixed_structure), size=16, replace=False)
    sequences_titin_fixed_structure = list(df_titin_fixed_structure['peptide'].values[indices])
    while len(set(sequences_titin_fixed_structure).intersection(set(all_sequences))) > 0 or len(set(sequences_titin_fixed_structure)) < len(sequences_titin_fixed_structure):
        indices = np.random.choice(len(df_titin_fixed_structure), size=16, replace=False)
        sequences_titin_fixed_structure = list(df_titin_fixed_structure['peptide'].values[indices])
    paes_titin_fixed_structure = list(df_titin_fixed_structure['pmhc_tcr_pae'].values[indices])
    all_sequences.extend(sequences_titin_fixed_structure)
    all_paes.extend(paes_titin_fixed_structure)
    all_model_types.extend(list(df_titin_fixed_structure['hcnn_model'].values[indices]))

    df_titin_sim_anneal = df_titin_sim_anneal.loc[np.logical_and(df_titin_sim_anneal['is_binder_by_netmhc_pan'] == 1, df_titin_sim_anneal['pmhc_tcr_pae'] < PAE_CUTOFF)]
    min_hamming_distances = np.array([min([hamming_distance(seq, WT_SEQ) for WT_SEQ in sequences_wildtype]) for seq in df_titin_sim_anneal['peptide'].values])
    df_titin_sim_anneal = df_titin_sim_anneal.loc[min_hamming_distances >= 3]
    indices = np.random.choice(len(df_titin_sim_anneal), size=16, replace=False)
    sequences_titin_sim_anneal = list(df_titin_sim_anneal['peptide'].values[indices])
    while len(set(sequences_titin_sim_anneal).intersection(set(all_sequences))) > 0 or len(set(sequences_titin_sim_anneal)) < len(sequences_titin_sim_anneal):
        indices = np.random.choice(len(df_titin_sim_anneal), size=16, replace=False)
        sequences_titin_sim_anneal = list(df_titin_sim_anneal['peptide'].values[indices])
    paes_titin_sim_anneal = list(df_titin_sim_anneal['pmhc_tcr_pae'].values[indices])
    all_sequences.extend(sequences_titin_sim_anneal)
    all_paes.extend(paes_titin_sim_anneal)
    all_model_types.extend(list(df_titin_sim_anneal['hcnn_model'].values[indices]))


    df_mage_with_e1_fixed_structure = df_mage_with_e1_fixed_structure.loc[np.logical_and(df_mage_with_e1_fixed_structure['is_binder_by_netmhc_pan'] == 1, df_mage_with_e1_fixed_structure['pmhc_tcr_pae'] < PAE_CUTOFF)]
    min_hamming_distances = np.array([min([hamming_distance(seq, WT_SEQ) for WT_SEQ in sequences_wildtype]) for seq in df_mage_with_e1_fixed_structure['peptide'].values])
    df_mage_with_e1_fixed_structure = df_mage_with_e1_fixed_structure.loc[min_hamming_distances >= 3]
    indices = np.random.choice(len(df_mage_with_e1_fixed_structure), size=15, replace=False)
    sequences_mage_with_e1_fixed_structure = list(df_mage_with_e1_fixed_structure['peptide'].values[indices])
    while len(set(sequences_mage_with_e1_fixed_structure).intersection(set(all_sequences))) > 0 or len(set(sequences_mage_with_e1_fixed_structure)) < len(sequences_mage_with_e1_fixed_structure):
        indices = np.random.choice(len(df_mage_with_e1_fixed_structure), size=15, replace=False)
        sequences_mage_with_e1_fixed_structure = list(df_mage_with_e1_fixed_structure['peptide'].values[indices])
    paes_mage_with_e1_fixed_structure = list(df_mage_with_e1_fixed_structure['pmhc_tcr_pae'].values[indices])
    all_sequences.extend(sequences_mage_with_e1_fixed_structure)
    all_paes.extend(paes_mage_with_e1_fixed_structure)
    all_model_types.extend([model + '__fixed_E1' for model in list(df_mage_with_e1_fixed_structure['hcnn_model'].values[indices])])


    df_mage_with_e1_sim_anneal = df_mage_with_e1_sim_anneal.loc[np.logical_and(df_mage_with_e1_sim_anneal['is_binder_by_netmhc_pan'] == 1, df_mage_with_e1_sim_anneal['pmhc_tcr_pae'] < PAE_CUTOFF)]
    min_hamming_distances = np.array([min([hamming_distance(seq, WT_SEQ) for WT_SEQ in sequences_wildtype]) for seq in df_mage_with_e1_sim_anneal['peptide'].values])
    df_mage_with_e1_sim_anneal = df_mage_with_e1_sim_anneal.loc[min_hamming_distances >= 3]
    indices = np.random.choice(len(df_mage_with_e1_sim_anneal), size=15, replace=False)
    sequences_mage_with_e1_sim_anneal = list(df_mage_with_e1_sim_anneal['peptide'].values[indices])
    while len(set(sequences_mage_with_e1_sim_anneal).intersection(set(all_sequences))) > 0 or len(set(sequences_mage_with_e1_sim_anneal)) < len(sequences_mage_with_e1_sim_anneal):
        indices = np.random.choice(len(df_mage_with_e1_sim_anneal), size=15, replace=False)
        sequences_mage_with_e1_sim_anneal = list(df_mage_with_e1_sim_anneal['peptide'].values[indices])
    paes_mage_with_e1_sim_anneal = list(df_mage_with_e1_sim_anneal['pmhc_tcr_pae'].values[indices])
    all_sequences.extend(sequences_mage_with_e1_sim_anneal)
    all_paes.extend(paes_mage_with_e1_sim_anneal)
    all_model_types.extend([model + '__fixed_E1' for model in list(df_mage_with_e1_sim_anneal['hcnn_model'].values[indices])])

    assert len(all_sequences) == len(set(all_sequences)), f'Duplicate sequences! {len(all_sequences)} vs. {len(set(all_sequences))}'

    sequence_sets = [sequences_wildtype, sequences_mage_fixed_structure, sequences_mage_sim_anneal, sequences_titin_fixed_structure, sequences_titin_sim_anneal, sequences_mage_with_e1_fixed_structure, sequences_mage_with_e1_sim_anneal]
    pae_sets = [paes_wildtype, paes_mage_fixed_structure, paes_mage_sim_anneal, paes_titin_fixed_structure, paes_titin_sim_anneal, paes_mage_with_e1_fixed_structure, paes_mage_with_e1_sim_anneal]
    

    ## plot hamming distances to WT
    plt.figure(figsize=(22, 5))

    tab20_colors = [cm.to_hex(plt.cm.tab20(i)) for i in range(20)]
    color_sets = ['black', tab20_colors[6], tab20_colors[7], tab20_colors[0], tab20_colors[1], tab20_colors[12], tab20_colors[13]]
    color_names = ['Wildtype', 'MageA3 Fixed Struc.', 'MageA3 Sim. Anneal.', 'Titin Fixed Struc.' 'Titin Sim. Anneal.', 'MageA3 Fixed Struc. with E1', 'MageA3 Sim. Anneal. with E1']

    plt.grid(axis='y', ls='--', alpha=0.5)
    
    start_ind = 0
    for sequences, color in zip(sequence_sets,
                                color_sets):
        indices = np.arange(start_ind, start_ind+len(sequences))
        start_ind += len(sequences)

        for WT_SEQ, marker in zip([MAGE_SEQ, TITIN_SEQ],
                                  ['s', '^']):
            hamming_distances = [hamming_distance(seq, WT_SEQ) for seq in sequences]
            plt.scatter(indices, hamming_distances, color=color, marker=marker)
    
    plt.xticks(np.arange(len(all_sequences)), ['MageA3', 'Titin'] + [f'p{i+1}' for i in range(len(all_sequences)-2)], rotation=80, fontsize=11)
    plt.yticks(np.arange(10), np.arange(10), fontsize=11)

    plt.ylabel('Hamming Distance\nfrom MageA3 and Titin', fontsize=14)

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    legend_handles = [Patch(edgecolor='w', facecolor=color, label=name) for color, name in zip(color_sets, color_names)]
    legend_handles += [Line2D([0], [0], marker='s', color='w', label='Distance from MageA3', markerfacecolor='k', markersize=15),
                        Line2D([0], [0], marker='^', color='w', label='Distance from Titin', markerfacecolor='k', markersize=15)]
    plt.gca().legend(handles=legend_handles, loc='center right', fontsize=11, bbox_to_anchor=(1.2, 0.5))

    plt.tight_layout()
    plt.savefig(f'magea3_and_titin_hamming_distance_of_selected_designs__pae_cutoff={PAE_CUTOFF}__seed={SEED}.png')
    plt.close()


    ## plot PAE
    plt.figure(figsize=(22, 5))

    tab20_colors = [cm.to_hex(plt.cm.tab20(i)) for i in range(20)]
    color_sets = ['black', tab20_colors[6], tab20_colors[7], tab20_colors[0], tab20_colors[1], tab20_colors[12], tab20_colors[13]]
    color_names = ['Wildtype', 'MageA3 Fixed Struc.', 'MageA3 Sim. Anneal.', 'Titin Fixed Struc.' 'Titin Sim. Anneal.', 'MageA3 Fixed Struc. with E1', 'MageA3 Sim. Anneal. with E1']

    plt.grid(axis='y', ls='--', alpha=0.5)
    
    start_ind = 0
    for sequences, paes, color in zip(sequence_sets,
                                      pae_sets,
                                      color_sets):
        indices = np.arange(start_ind, start_ind+len(sequences))
        start_ind += len(sequences)

        plt.scatter(indices, -np.array(paes), color=color, marker='o')
    
    plt.xticks(np.arange(len(all_sequences)), ['MageA3', 'Titin'] + [f'p{i+1}' for i in range(len(all_sequences)-2)], rotation=80, fontsize=11)

    plt.ylabel('Negative TCRDock PAE', fontsize=14)

    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch

    legend_handles = [Patch(edgecolor='w', facecolor=color, label=name) for color, name in zip(color_sets, color_names)]
    plt.gca().legend(handles=legend_handles, loc='center right', fontsize=11, bbox_to_anchor=(1.2, 0.5))

    plt.tight_layout()
    plt.savefig(f'magea3_and_titin_neg_tcrdock_pae_of_selected_designs__pae_cutoff={PAE_CUTOFF}__seed={SEED}.png')
    plt.close()


    ## save peptides
    wt_for_design = []
    wt_for_design += ['magea3', 'titin']
    wt_for_design += ['magea3'] * len(sequences_mage_fixed_structure)
    wt_for_design += ['magea3'] * len(sequences_mage_sim_anneal)
    wt_for_design += ['titin'] * len(sequences_titin_fixed_structure)
    wt_for_design += ['titin'] * len(sequences_titin_sim_anneal)
    wt_for_design += ['magea3'] * len(sequences_mage_with_e1_fixed_structure)
    wt_for_design += ['magea3'] * len(sequences_mage_with_e1_sim_anneal)

    sampling_method = []
    sampling_method += ['N/A', 'N/A']
    sampling_method += ['fixed_structure'] * len(sequences_mage_fixed_structure)
    sampling_method += ['sim_anneal'] * len(sequences_mage_sim_anneal)
    sampling_method += ['fixed_structure'] * len(sequences_titin_fixed_structure)
    sampling_method += ['sim_anneal'] * len(sequences_titin_sim_anneal)
    sampling_method += ['fixed_structure'] * len(sequences_mage_with_e1_fixed_structure)
    sampling_method += ['sim_anneal'] * len(sequences_mage_with_e1_sim_anneal)

    hamming_distances_from_magea3 = [hamming_distance(seq, MAGE_SEQ) for seq in all_sequences]
    hamming_distances_from_titin = [hamming_distance(seq, TITIN_SEQ) for seq in all_sequences]

    df = pd.DataFrame({
        'sequence': all_sequences,
        'model': all_model_types,
        'sampling_method': sampling_method,
        'wt_for_design': wt_for_design,
        'tcrdock_pae': all_paes,
        'hamming_distances_from_magea3': hamming_distances_from_magea3,
        'hamming_distance_from_titin': hamming_distances_from_titin
    })

    df.to_csv(f'magea3_and_titin_peptide_designs_with_hcnn__pae_cutoff={PAE_CUTOFF}__seed={SEED}.csv', index=None)


    






