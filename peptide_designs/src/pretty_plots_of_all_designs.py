

import os, sys
import numpy as np
import pandas as pd

from scipy.stats import pearsonr, spearmanr

import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import argparse

THIS_FILE = os.path.dirname(os.path.abspath(__file__))

def get_dfs_for_system(system):

    if system == 'nyeso':

        df_wildtype = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/wildtype/wildtype_w_pae_w_blosum.tsv'), sep='\t')

        df_hermes_relaxed_000 = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv'), sep='\t')
        df_hermes_relaxed_050 = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv'), sep='\t')

        df_mhc = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/mhc_pwm/mhc_motif_peptides_w_pae_w_blosum.tsv'), sep='\t')

        df_proteinmpnn_002 = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7_w_pae_w_blosum.tsv'), sep='\t')
        df_proteinmpnn_020 = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7_w_pae_w_blosum.tsv'), sep='\t')

        df_hermes_fixed_000 = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv'), sep='\t')
        df_hermes_fixed_050 = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv'), sep='\t')
        
        df_blosum_t1 = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/blosum62/sample_peptides_from_blosum62__temperature=1.0_w_pae_w_blosum.tsv'), sep='\t')
        df_blosum_t2 = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae_w_blosum.tsv'), sep='\t')
        df_blosum_t3 = pd.read_csv(os.path.join(THIS_FILE, '../nyeso/nyeso_full_copy/blosum62/sample_peptides_from_blosum62__temperature=3.0_w_pae_w_blosum.tsv'), sep='\t')

    elif system == 'magea3_and_titin':

        df_wildtype = None
        df_hermes_relaxed_000 = None
        df_hermes_relaxed_050 = None
        df_mhc = None
        df_proteinmpnn_002 = None
        df_proteinmpnn_020 = None
        df_hermes_fixed_000 = None
        df_hermes_fixed_050 = None
        df_blosum_t1 = None
        df_blosum_t2 = None
        df_blosum_t3 = None

        for i, substruct in enumerate(['magea3', 'magea3_with_e1', 'titin']):

            if df_wildtype is None or i < 2: # hack to keep only one copy of magea3 sequence
                df_wildtype = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/wildtype/wildtype_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_wildtype = pd.concat([df_wildtype, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/wildtype/wildtype_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_hermes_relaxed_000 is None:
                df_hermes_relaxed_000 = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_hermes_relaxed_000 = pd.concat([df_hermes_relaxed_000, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_hermes_relaxed_050 is None:
                df_hermes_relaxed_050 = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_hermes_relaxed_050 = pd.concat([df_hermes_relaxed_050, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_mhc is None:
                df_mhc = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/mhc_pwm/mhc_motif_peptides_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_mhc = pd.concat([df_mhc, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/mhc_pwm/mhc_motif_peptides_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_proteinmpnn_002 is None:
                df_proteinmpnn_002 = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_proteinmpnn_002 = pd.concat([df_proteinmpnn_002, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7_w_pae_w_blosum.tsv'), sep='\t')])

            if df_proteinmpnn_020 is None:
                df_proteinmpnn_020 = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_proteinmpnn_020 = pd.concat([df_proteinmpnn_020, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_hermes_fixed_000 is None:
                df_hermes_fixed_000 = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_hermes_fixed_000 = pd.concat([df_hermes_fixed_000, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_hermes_fixed_050 is None:
                df_hermes_fixed_050 = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_hermes_fixed_050 = pd.concat([df_hermes_fixed_050, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_blosum_t1 is None:
                df_blosum_t1 = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=1.0_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_blosum_t1 = pd.concat([df_blosum_t1, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=1.0_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_blosum_t2 is None:
                df_blosum_t2 = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_blosum_t2 = pd.concat([df_blosum_t2, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_blosum_t3 is None:
                df_blosum_t3 = pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=3.0_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_blosum_t3 = pd.concat([df_blosum_t3, pd.read_csv(os.path.join(THIS_FILE, f'../magea3_and_titin/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=3.0_w_pae_w_blosum.tsv'), sep='\t')])

    elif system == 'ebv':

        df_wildtype = None
        df_hermes_relaxed_000 = None
        df_hermes_relaxed_050 = None
        df_mhc = None
        df_proteinmpnn_002 = None
        df_proteinmpnn_020 = None
        df_hermes_fixed_000 = None
        df_hermes_fixed_050 = None
        df_blosum_t1 = None
        df_blosum_t2 = None
        df_blosum_t3 = None

        for substruct in ['hpvg', 'hpvg_q5']:

            if df_wildtype is None:
                df_wildtype = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/wildtype/wildtype_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_wildtype = pd.concat([df_wildtype, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/wildtype/wildtype_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_hermes_relaxed_000 is None:
                df_hermes_relaxed_000 = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_hermes_relaxed_000 = pd.concat([df_hermes_relaxed_000, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_hermes_relaxed_050 is None:
                df_hermes_relaxed_050 = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_hermes_relaxed_050 = pd.concat([df_hermes_relaxed_050, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_mhc is None:
                df_mhc = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/mhc_pwm/mhc_motif_peptides_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_mhc = pd.concat([df_mhc, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/mhc_pwm/mhc_motif_peptides_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_proteinmpnn_002 is None:
                df_proteinmpnn_002 = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_proteinmpnn_002 = pd.concat([df_proteinmpnn_002, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7_w_pae_w_blosum.tsv'), sep='\t')])

            if df_proteinmpnn_020 is None:
                df_proteinmpnn_020 = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_proteinmpnn_020 = pd.concat([df_proteinmpnn_020, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_hermes_fixed_000 is None:
                df_hermes_fixed_000 = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_hermes_fixed_000 = pd.concat([df_hermes_fixed_000, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_hermes_fixed_050 is None:
                df_hermes_fixed_050 = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_hermes_fixed_050 = pd.concat([df_hermes_fixed_050, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_blosum_t1 is None:
                df_blosum_t1 = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=1.0_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_blosum_t1 = pd.concat([df_blosum_t1, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=1.0_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_blosum_t2 is None:
                df_blosum_t2 = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_blosum_t2 = pd.concat([df_blosum_t2, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae_w_blosum.tsv'), sep='\t')])
            
            if df_blosum_t3 is None:
                df_blosum_t3 = pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=3.0_w_pae_w_blosum.tsv'), sep='\t')
            else:
                df_blosum_t3 = pd.concat([df_blosum_t3, pd.read_csv(os.path.join(THIS_FILE, f'../ebv/{substruct}/blosum62/sample_peptides_from_blosum62__temperature=3.0_w_pae_w_blosum.tsv'), sep='\t')])

    # add min hamming distance as extra column
    # dfs = [df_hermes_relaxed_000, df_hermes_relaxed_050, df_mhc, df_proteinmpnn_002, df_proteinmpnn_020, df_hermes_fixed_000, df_hermes_fixed_050, df_blosum_t1, df_blosum_t2, df_blosum_t3]
    wildtypes = df_wildtype['peptide'].values

    def hamming_distance(seq1, seq2):
        assert len(seq1) == len(seq2)
        return sum([1 for i in range(len(seq1)) if seq1[i] != seq2[i]])
    
    def add_hamming_distances_to_df(df, wildtypes):
        wt_to_hamming = {f'hamming_distance_to_{wt_seq}': [] for wt_seq in wildtypes}
        wt_to_hamming['min_hamming_distance_to_wt'] = []
        for i, row in df.iterrows():
            curr_seq = row['peptide']
            min_val = np.inf
            for wt_seq in wildtypes:
                curr_dist = hamming_distance(curr_seq, wt_seq)
                wt_to_hamming[f'hamming_distance_to_{wt_seq}'].append(curr_dist)
                min_val = min(min_val, curr_dist)
            wt_to_hamming['min_hamming_distance_to_wt'].append(min_val)
        for entry in wt_to_hamming:
            df[entry] = wt_to_hamming[entry]
        return df

    def remove_duplicates(df): # I forgot to do it for blosum designs, and also there may be duplicates when merging designs from different base structure
        return df.drop_duplicates(subset=['peptide'], keep='first', inplace=False)

    df_hermes_relaxed_000 = remove_duplicates(add_hamming_distances_to_df(df_hermes_relaxed_000, wildtypes))
    df_hermes_relaxed_050 = remove_duplicates(add_hamming_distances_to_df(df_hermes_relaxed_050, wildtypes))
    df_mhc = remove_duplicates(add_hamming_distances_to_df(df_mhc, wildtypes))
    df_proteinmpnn_002 = remove_duplicates(add_hamming_distances_to_df(df_proteinmpnn_002, wildtypes))
    df_proteinmpnn_020 = remove_duplicates(add_hamming_distances_to_df(df_proteinmpnn_020, wildtypes))
    df_hermes_fixed_000 = remove_duplicates(add_hamming_distances_to_df(df_hermes_fixed_000, wildtypes))
    df_hermes_fixed_050 = remove_duplicates(add_hamming_distances_to_df(df_hermes_fixed_050, wildtypes))
    df_blosum_t1 = remove_duplicates(add_hamming_distances_to_df(df_blosum_t1, wildtypes))
    df_blosum_t2 = remove_duplicates(add_hamming_distances_to_df(df_blosum_t2, wildtypes))
    df_blosum_t3 = remove_duplicates(add_hamming_distances_to_df(df_blosum_t3, wildtypes))

    return {'wt': df_wildtype, 'hermes_relaxed_000': df_hermes_relaxed_000, 'hermes_relaxed_050': df_hermes_relaxed_050, 'mhc': df_mhc, 'proteinmpnn_002': df_proteinmpnn_002, 'proteinmpnn_020': df_proteinmpnn_020, 'hermes_fixed_000': df_hermes_fixed_000, 'hermes_fixed_050': df_hermes_fixed_050, 'blosum_t1': df_blosum_t1, 'blosum_t2': df_blosum_t2, 'blosum_t3': df_blosum_t3}

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
    'hermes_relaxed_000': orange,
    'hermes_relaxed_050': orange_light,
    'hermes_fixed_000': red,
    'hermes_fixed_050': red_light,
    'mhc': gray,
    'proteinmpnn_002': blue,
    'proteinmpnn_020': blue_light,
    'blosum_t1': brown,
    'blosum_t2': brown_light,
    'blosum_t3': '#dec6c2',
}

MODEL_TO_PRETTY_NAME = {
    'hermes_relaxed_000': 'HERMES-$relaxed$ 0.00',
    'hermes_relaxed_050': 'HERMES-$relaxed$ 0.50',
    'hermes_fixed_000': 'HERMES-$fixed$ 0.00',
    'hermes_fixed_050': 'HERMES-$fixed$ 0.50',
    'mhc': '{mhc} PWM',
    'proteinmpnn_002': 'ProteinMPNN 0.02',
    'proteinmpnn_020': 'ProteinMPNN 0.20',
    'blosum_t1': 'BLOSUM62 temp = 1.0',
    'blosum_t2': 'BLOSUM62 temp = 2.0',
    'blosum_t3': 'BLOSUM62 temp = 3.0'
}

SYSTEM_TO_MHC = {
    'nyeso': 'A*02:01',
    'magea3_and_titin': 'A*01:01',
    'ebv': 'B*35:01'
}

SYSTEM_TO_WILDTYPE_NAMES = {
    'nyeso': ['NYESO', 'NYESO C9V'],
    'magea3_and_titin': ['MageA3', 'Titin'],
    'ebv': ['HPVG', 'HPVG E5Q']
}

SYSTEM_TO_PEP_LENGTH = {
    'nyeso': 9,
    'magea3_and_titin': 9,
    'ebv': 11
}

MODELS_IN_ORDER = ['hermes_fixed_000', 'hermes_fixed_050', 'hermes_relaxed_000', 'hermes_relaxed_050', 'proteinmpnn_002', 'proteinmpnn_020', 'blosum_t1', 'blosum_t2', 'blosum_t3', 'mhc']

SYSTEM_TO_PAE_THRESHOLD = {
    'nyeso': 5.5,
    'magea3_and_titin': 5.2,
    'ebv': 5.1
}

ALPHA = 0.5

def color_violinplot_multiple(parts, colors):
    for pm, pb, color in zip(parts['cmedians'], parts['bodies'], colors):
        pm.set_color(color)
        pb.set_facecolor(color)
        pb.set_edgecolor(color)
        pb.set_alpha(ALPHA)

def color_violinplot_single(parts, color):
    parts['cmedians'].set_color(color)
    parts['cmedians'].set_linewidth(2)
    
    for pb, color in zip(parts['bodies'], [color]):
        pb.set_facecolor(color)
        pb.set_edgecolor(color)
        pb.set_alpha(ALPHA)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--system', type=str, required=True)
    args = parser.parse_args()


    model_to_df = get_dfs_for_system(args.system)


    ## PAE, hamming distance, is_mhc_binder

    fontsize = 20

    ncols = 4
    nrows = 1
    colsize = 5
    rowsize = 5
    fig, axs = plt.subplots(figsize=(ncols*colsize, nrows*rowsize), ncols=ncols, nrows=nrows)


    ## 1) negative PAE violinplot/boxplot

    ax = axs[0]

    positions = range(len(MODELS_IN_ORDER))

    for pos, model in zip(positions, MODELS_IN_ORDER):
        
        scores = model_to_df[model]['neg_pmhc_tcr_pae'].values
        color = MODEL_TO_COLOR[model]

        parts = ax.violinplot([scores], positions=[pos], widths=0.6, showmedians=True, showextrema=False)
        color_violinplot_single(parts, [color])
    
    # put number of examples under each violinplot
    # do it in separate loop so that the ylims get updated to their final values
    # UPDATE: showing this in a separate plot
    # for pos, model in zip(positions, MODELS_IN_ORDER):
        
    #     scores = model_to_df[model]['neg_pmhc_tcr_pae'].values
    #     num_scores = len(scores)

    #     min_ylim, max_ylim = ax.get_ylim()
    #     ax.text(pos, min_ylim-(0.06 + 0.04*((pos%2)-1))*min_ylim, num_scores, ha='center', va='bottom', fontsize=fontsize-5)
    
    # wildtype values, use square and diamond like in the main text figure, the order should be preserved, keep it simple without dictionaries
    markers = ['s', 'd']
    wt_pos = min(positions) - 1
    for i, (marker, wt_score) in enumerate(zip(markers, model_to_df['wt']['neg_pmhc_tcr_pae'])):
        ax.scatter(wt_pos+(i-1)*0.3, wt_score, color='black', marker=marker)

    
    # horizontal line with our decided pae threshold
    ax.axhline(-SYSTEM_TO_PAE_THRESHOLD[args.system], ls='--', color='black')

    ax.set_xticks([wt_pos], ['WT'])
    ax.tick_params(axis='both', labelsize=fontsize-2)
    ax.set_ylabel('- TCRdock PAE', fontsize=fontsize)



    ## 2) min hamming distance heatmap

    ax = axs[1]

    positions = range(len(MODELS_IN_ORDER))

    # for pos, model in zip(positions, MODELS_IN_ORDER):
        
    #     hamming_distances = model_to_df[model]['min_hamming_distance_to_wt'].values
    #     color = MODEL_TO_COLOR[model]

    #     parts = ax.violinplot([hamming_distances], positions=[pos], widths=0.6, showmedians=True, showextrema=False)
    #     color_violinplot_single(parts, [color])
    
    # ax.grid(axis='y', ls='--', color='dimgrey', alpha=0.5)

    # ax.set_xticks([])
    # ax.tick_params(axis='both', labelsize=fontsize-2)
    # ax.set_ylabel('hamming dist. from WT', fontsize=fontsize)

    import matplotlib.colors as mcolors
    from collections import Counter

    def custom_colormap(color):
        return mcolors.LinearSegmentedColormap.from_list("custom_cmap", ["white", color])
    
    all_hamming_distances = np.arange(SYSTEM_TO_PEP_LENGTH[args.system]+1)

    
    x = np.arange(len(MODELS_IN_ORDER)+1)
    y = np.arange(len(all_hamming_distances)+1)

    for i_col, model in zip(positions, MODELS_IN_ORDER):
        
        hamming_distances = model_to_df[model]['min_hamming_distance_to_wt'].values
        color = MODEL_TO_COLOR[model]

        counts = dict(Counter(hamming_distances))


        h_dist_counts_column = []
        for h_dist in all_hamming_distances:
            if h_dist not in counts:
                h_dist_counts_column.append(0)
            else:
                h_dist_counts_column.append(counts[h_dist])
        
        h_dist_counts_column = np.array(h_dist_counts_column)

        # Create a mesh grid for this column (taking one slice of x at a time)
        X, Y = np.meshgrid(x[i_col:i_col+2], y)

        # Normalize data for each column
        cmap_norm = mcolors.Normalize(vmin=np.min(h_dist_counts_column), vmax=np.max(h_dist_counts_column))
        
        # Plot the column using the correct shape
        ax.pcolormesh(X, Y, h_dist_counts_column.reshape(-1, 1), cmap=custom_colormap(color), norm=cmap_norm)

        # add number in each cell
        for i_row_and_h_dist in all_hamming_distances:
            if h_dist_counts_column[i_row_and_h_dist] != 0:
                ax.text(i_col+0.5, i_row_and_h_dist+0.5, h_dist_counts_column[i_row_and_h_dist], ha='center', va='center', color='black', fontsize=13)
        
        # add total number below the column
        ax.text(i_col+0.5, -0.5 - 0.5*(i_col % 2), int(np.sum(h_dist_counts_column)), ha='center', va='center', color='black', fontsize=13)
    
    # ax.text(-0.25, -0.75, 'count', ha='right', va='center', color='black', fontsize=fontsize-3)
    
    ax.set_ylabel('hamming dist. from WT', fontsize=fontsize)
    ax.set_xticks([])
    ax.set_yticks(all_hamming_distances+0.5, all_hamming_distances, fontsize=fontsize-2)
    ax.set_title('number of designs', fontsize=fontsize)


    ## 3) min hamming distance heatmap with fraction of designs

    ax = axs[2]

    all_hamming_distances = np.arange(SYSTEM_TO_PEP_LENGTH[args.system]+1)

    
    x = np.arange(len(MODELS_IN_ORDER)+1)
    y = np.arange(len(all_hamming_distances)+1)

    for i_col, model in zip(positions, MODELS_IN_ORDER):
        
        hamming_distances = model_to_df[model]['min_hamming_distance_to_wt'].values
        is_above_threshold = model_to_df[model]['neg_pmhc_tcr_pae'].values > -SYSTEM_TO_PAE_THRESHOLD[args.system]
        color = MODEL_TO_COLOR[model]


        h_dist_counts_column = []
        for h_dist in all_hamming_distances:
            if h_dist not in hamming_distances:
                h_dist_counts_column.append(np.nan)
            else:
                h_dist_counts_column.append(int(np.mean(is_above_threshold[hamming_distances == h_dist])*100))
        
        h_dist_counts_column = np.array(h_dist_counts_column)

        # Create a mesh grid for this column (taking one slice of x at a time)
        X, Y = np.meshgrid(x[i_col:i_col+2], y)

        # Normalize data globally
        cmap_norm = mcolors.Normalize(vmin=0, vmax=100)
        
        # Plot the column using the correct shape
        ax.pcolormesh(X, Y, h_dist_counts_column.reshape(-1, 1), cmap=custom_colormap(color), norm=cmap_norm)

        # add number in each cell
        for i_row_and_h_dist in all_hamming_distances:
            if not np.isnan(h_dist_counts_column[i_row_and_h_dist]):
                ax.text(i_col+0.5, i_row_and_h_dist+0.5, int(h_dist_counts_column[i_row_and_h_dist]), ha='center', va='center', color='black', fontsize=13)
        
        # # add total number below the column
        # ax.text(i_col+0.5, -0.5 - 0.5*(i_col % 2), int(np.mean(is_above_threshold)*100), ha='center', va='center', color='black', fontsize=13)
    
    # ax.text(-0.25, -0.75, 'count', ha='right', va='center', color='black', fontsize=fontsize-3)
    
    ax.set_ylabel('hamming dist. from WT', fontsize=fontsize)
    ax.set_xticks([])
    ax.set_yticks(all_hamming_distances+0.5, all_hamming_distances, fontsize=fontsize-2)
    ax.set_title('fraction above\n- TCRdock PAE threshold (%)', fontsize=fontsize)



    ## 4) barplot of is_mhc_binder

    ax = axs[3]

    positions = range(len(MODELS_IN_ORDER))

    for pos, model in zip(positions, MODELS_IN_ORDER):
        
        is_mhc_binder_weak = model_to_df[model]['is_binder_by_netmhcpan__weak'].values
        is_mhc_binder_strong = model_to_df[model]['is_binder_by_netmhcpan'].values

        perc_weak = np.mean(is_mhc_binder_weak) * 100
        perc_strong = np.mean(is_mhc_binder_strong) * 100
        perc_weak_but_not_strong = perc_weak - perc_strong

        color = MODEL_TO_COLOR[model]

        ax.bar([pos], [perc_strong], color=color, hatch='...', alpha=ALPHA)
        ax.bar([pos], [perc_weak_but_not_strong], bottom=[perc_strong], color=color, alpha=ALPHA)
    
    ax.grid(axis='y', ls='--', alpha=0.5)
    
    ax.set_ylim([0, 100])

    ax.set_xticks([])
    ax.tick_params(axis='both', labelsize=fontsize-2)
    ax.set_ylabel(f'fraction {SYSTEM_TO_MHC[args.system]} binder (%)', fontsize=fontsize)

    plt.tight_layout()
    plt.savefig(f'../{args.system}/design_summary_plot.png')
    plt.savefig(f'../{args.system}/design_summary_plot.pdf')
    plt.close()


    ## make legend
    handles = []
    for model in MODELS_IN_ORDER:
        color = MODEL_TO_COLOR[model]
        model_pretty_name = MODEL_TO_PRETTY_NAME[model]
        if model == 'mhc': model_pretty_name = model_pretty_name.format(mhc=SYSTEM_TO_MHC[args.system])
        handles.append(Patch(facecolor=color, edgecolor=color, alpha=ALPHA, label=model_pretty_name))
    for marker, wt_name in zip(['s', 'd'], SYSTEM_TO_WILDTYPE_NAMES[args.system]):
        handles.append(Line2D([0], [0], color='black', linestyle='', marker=marker, label=wt_name))
    handles.append(Line2D([0], [0], color='black', linestyle='--', label='chosen PAE threshold'))

    handles.append(Patch(facecolor='none', edgecolor='black', hatch='...', label='strong binder'))
    handles.append(Patch(facecolor='none', edgecolor='black', label='weak binder'))

    plt.figure(figsize=(5, 5))
    plt.legend(handles=handles, loc='center', fontsize=14)
    plt.axis('off')
    plt.savefig(f'../{args.system}/design_summary_plot__legend.png')
    plt.savefig(f'../{args.system}/design_summary_plot__legend.pdf')
    plt.close()

    plt.figure(figsize=(15, 2))
    plt.legend(handles=handles, loc='center', fontsize=14, ncol=len(handles)//3 + len(handles)%3)
    plt.axis('off')
    plt.savefig(f'../{args.system}/design_summary_plot__legend_horizontal.png')
    plt.savefig(f'../{args.system}/design_summary_plot__legend_horizontal.pdf')
    plt.close()



    ## plot of pae by hamming distance from wildtype

    min_hamming_distance = np.inf
    max_hamming_distance = -np.inf
    for model in MODELS_IN_ORDER:
        df = model_to_df[model]
        distances = df['min_hamming_distance_to_wt'].values
        min_hamming_distance = min(min_hamming_distance, min(distances))
        max_hamming_distance = max(min_hamming_distance, max(distances))

    hamming_distances = range(min_hamming_distance, max_hamming_distance+1)

    nrows = 2
    ncols = (len(hamming_distances) // 2) + (len(hamming_distances) % 2)
    colsize = 3.5
    rowsize = 3.5
    fig, axs = plt.subplots(figsize=(ncols*colsize, nrows*rowsize), ncols=ncols, nrows=nrows, sharex=True, sharey=True)


    for i, h_dist in enumerate(hamming_distances):

        ax = axs.flatten()[i]

        positions = range(len(MODELS_IN_ORDER))

        for pos, model in zip(positions, MODELS_IN_ORDER):
            
            scores = model_to_df[model]['neg_pmhc_tcr_pae'].values
            h_dist_mask = model_to_df[model]['min_hamming_distance_to_wt'].values == h_dist
            scores = scores[h_dist_mask]

            if len(scores) == 0: continue

            color = MODEL_TO_COLOR[model]

            parts = ax.violinplot([scores], positions=[pos], widths=0.6, showmedians=True, showextrema=False)
            color_violinplot_single(parts, [color])

        ax.set_title(f'hamming dist. = {h_dist}', fontsize=fontsize-2)
        ax.set_xticks([])

        # horizontal line with our decided pae threshold
        ax.axhline(-SYSTEM_TO_PAE_THRESHOLD[args.system], ls='--', color='black')
    
    # put number of examples under each violinplot
    # do it in separate loop so that the ylims get updated to their final values
    for i, h_dist in enumerate(hamming_distances):

        ax = axs.flatten()[i]

        positions = range(len(MODELS_IN_ORDER))
    
        for pos, model in zip(positions, MODELS_IN_ORDER):
            
            scores = model_to_df[model]['neg_pmhc_tcr_pae'].values
            h_dist_mask = model_to_df[model]['min_hamming_distance_to_wt'].values == h_dist
            scores = scores[h_dist_mask]

            if len(scores) == 0: continue

            num_scores = len(scores)
            min_ylim, max_ylim = ax.get_ylim()
            ax.text(pos, min_ylim-(0.06 + 0.04*((pos%2)-1))*min_ylim, num_scores, ha='center', va='bottom', fontsize=fontsize-6)
            
    for col in [0, ncols]:
        axs.flatten()[col].tick_params(axis='both', labelsize=fontsize-3)
        axs.flatten()[col].set_ylabel('- TCRdock PAE', fontsize=fontsize)

    plt.tight_layout()
    plt.savefig(f'../{args.system}/design_summary_pae_by_hamming_distance_from_wt.png')
    plt.savefig(f'../{args.system}/design_summary_pae_by_hamming_distance_from_wt.pdf')
    plt.close()

    
    fontsize = 22


    plt.figure(figsize=(8, 5))
    for model in MODELS_IN_ORDER:

        median_score_per_hamming_dist = []
        color = MODEL_TO_COLOR[model]

        for i, h_dist in enumerate(hamming_distances):
            
            scores = model_to_df[model]['neg_pmhc_tcr_pae'].values
            h_dist_mask = model_to_df[model]['min_hamming_distance_to_wt'].values == h_dist
            scores = scores[h_dist_mask]

            if len(scores) > 0:
                median_score_per_hamming_dist.append(np.median(scores))
            else:
                median_score_per_hamming_dist.append(np.nan)
    
        plt.plot(median_score_per_hamming_dist, marker='o', color=color, lw=3)

    plt.axhline(-SYSTEM_TO_PAE_THRESHOLD[args.system], ls='--', color='black')

    plt.xticks(hamming_distances, fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-3)
    plt.xlabel('hamming dist. from WT', fontsize=fontsize)
    plt.ylabel('- TCRdock PAE (median)', fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(f'../{args.system}/design_median_pae_vs_hamming_distance_from_wt.png')
    plt.savefig(f'../{args.system}/design_median_pae_vs_hamming_distance_from_wt.pdf')
    plt.close()


    ## barplot of fraction of designs above chosen threshold

    plt.figure(figsize=(6.2, 5))
    fraction_above_threshold = []
    colors = []
    for model in MODELS_IN_ORDER:
        scores = model_to_df[model]['neg_pmhc_tcr_pae'].values
        is_above_threshold = scores > -SYSTEM_TO_PAE_THRESHOLD[args.system]
        fraction_above_threshold.append(np.mean(is_above_threshold)*100)
        colors.append(MODEL_TO_COLOR[model])
    
    plt.bar(np.arange(len(MODELS_IN_ORDER)), fraction_above_threshold, color=colors, alpha=ALPHA)
    plt.xticks([])
    plt.yticks(fontsize=fontsize-1)
    plt.ylim([0, 100])

    # # write total number of designs at the bottom of the barplot
    # ax = plt.gca()
    # for pos, model in enumerate(MODELS_IN_ORDER):
    #     scores = model_to_df[model]['neg_pmhc_tcr_pae'].values
    #     num_scores = len(scores)
    #     min_ylim, max_ylim = ax.get_ylim()
    #     ax.text(pos, min_ylim-(0.06 + 0.04*((pos%2)-1))*min_ylim, num_scores, ha='center', va='bottom', fontsize=fontsize-6)
    
    plt.grid(axis='y', ls='--', alpha=0.5)

    plt.ylabel('fraction above\n- TCRdock PAE threshold', fontsize=fontsize+1)
    plt.tight_layout()
    plt.savefig(f'../{args.system}/design_above_threshold.png')
    plt.savefig(f'../{args.system}/design_above_threshold.pdf')
    plt.close()


    ## plot of fraction of designs above chosen threshold, as a function of hamming distance from WT

    plt.figure(figsize=(8, 5))
    for model in MODELS_IN_ORDER:

        fraction_above_threshold = []
        color = MODEL_TO_COLOR[model]

        for i, h_dist in enumerate(hamming_distances):
            
            scores = model_to_df[model]['neg_pmhc_tcr_pae'].values
            is_above_threshold = scores > -SYSTEM_TO_PAE_THRESHOLD[args.system]
            h_dist_mask = model_to_df[model]['min_hamming_distance_to_wt'].values == h_dist
            is_above_threshold = is_above_threshold[h_dist_mask]

            if len(is_above_threshold) > 0:
                fraction_above_threshold.append(np.mean(is_above_threshold)*100)
            else:
                fraction_above_threshold.append(np.nan)
    
        plt.plot(fraction_above_threshold, marker='o', color=color, lw=3)
    
    plt.grid(axis='y', ls='--', alpha=0.5)

    plt.xticks(hamming_distances, fontsize=fontsize-2)
    plt.yticks(fontsize=fontsize-3)
    plt.ylim([-5, 105])
    plt.xlabel('hamming dist. from WT', fontsize=fontsize)
    plt.ylabel('fraction above\n- TCRdock PAE threshold (%)', fontsize=fontsize)
    plt.tight_layout()
    plt.savefig(f'../{args.system}/design_above_threshold_vs_hamming_distance_from_wt.png')
    plt.savefig(f'../{args.system}/design_above_threshold_vs_hamming_distance_from_wt.pdf')
    plt.close()


