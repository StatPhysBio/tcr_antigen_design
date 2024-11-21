
import os
import numpy as np
import pandas as pd

def hamming_distance(seq1, seq2):
    assert len(seq1) == len(seq2)
    return sum([seq1[i] != seq2[i] for i in range(len(seq1))])

HPVG_SEQ = 'HPVGEADYFEY'
HPVG_Q5_SEQ = 'HPVGQADYFEY'

SEED = 43
np.random.seed(SEED)

PAE_CUTOFF = 5.1

if __name__ == '__main__':

    ## load the peptides sampled by the HCNN models

    df_hpvg_wt = pd.read_csv('hpvg/wildtype/wildtype_w_pae_w_blosum.tsv', sep='\t')
    df_hpvg_q5_wt = pd.read_csv('hpvg_q5/wildtype/wildtype_w_pae_w_blosum.tsv', sep='\t')

    sequences_wildtype = [HPVG_SEQ, HPVG_Q5_SEQ]

    df_hpvg_fixed_structure = pd.concat([pd.read_csv('hpvg/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                            pd.read_csv('hpvg/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])
    df_hpvg_fixed_structure = df_hpvg_fixed_structure.loc[np.logical_and(df_hpvg_fixed_structure['is_binder_by_netmhc_pan'] == 1, df_hpvg_fixed_structure['pmhc_tcr_pae'] < PAE_CUTOFF)]
    min_hamming_distances = np.array([min([hamming_distance(seq, WT_SEQ) for WT_SEQ in sequences_wildtype]) for seq in df_hpvg_fixed_structure['peptide'].values])
    df_hpvg_fixed_structure = df_hpvg_fixed_structure.loc[min_hamming_distances >= 3]


    df_hpvg_sim_anneal = pd.concat([pd.read_csv('hpvg/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                            pd.read_csv('hpvg/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])
    df_hpvg_sim_anneal = df_hpvg_sim_anneal.loc[np.logical_and(df_hpvg_sim_anneal['is_binder_by_netmhc_pan'] == 1, df_hpvg_sim_anneal['pmhc_tcr_pae'] < PAE_CUTOFF)]
    min_hamming_distances = np.array([min([hamming_distance(seq, WT_SEQ) for WT_SEQ in sequences_wildtype]) for seq in df_hpvg_sim_anneal['peptide'].values])
    df_hpvg_sim_anneal = df_hpvg_sim_anneal.loc[min_hamming_distances >= 3]

    df_hpvg_q5_fixed_structure = pd.concat([pd.read_csv('hpvg_q5/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                        pd.read_csv('hpvg_q5/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])
    df_hpvg_q5_fixed_structure = df_hpvg_q5_fixed_structure.loc[np.logical_and(df_hpvg_q5_fixed_structure['is_binder_by_netmhc_pan'] == 1, df_hpvg_q5_fixed_structure['pmhc_tcr_pae'] < PAE_CUTOFF)]
    min_hamming_distances = np.array([min([hamming_distance(seq, WT_SEQ) for WT_SEQ in sequences_wildtype]) for seq in df_hpvg_q5_fixed_structure['peptide'].values])
    df_hpvg_q5_fixed_structure = df_hpvg_q5_fixed_structure.loc[min_hamming_distances >= 3]

    df_hpvg_q5_sim_anneal = pd.concat([pd.read_csv('hpvg_q5/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t'),
                                        pd.read_csv('hpvg_q5/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')])
    df_hpvg_q5_sim_anneal = df_hpvg_q5_sim_anneal.loc[np.logical_and(df_hpvg_q5_sim_anneal['is_binder_by_netmhc_pan'] == 1, df_hpvg_q5_sim_anneal['pmhc_tcr_pae'] < PAE_CUTOFF)]
    min_hamming_distances = np.array([min([hamming_distance(seq, WT_SEQ) for WT_SEQ in sequences_wildtype]) for seq in df_hpvg_q5_sim_anneal['peptide'].values])
    df_hpvg_q5_sim_anneal = df_hpvg_q5_sim_anneal.loc[min_hamming_distances >= 3]


    ## load the selevted and filtered peptides
    df_chosen_first = pd.read_csv('hpvg_and_hpvg_q5_peptide_designs_with_hcnn__pae_cutoff=5.1__seed=43.csv')
    df_filtered = pd.read_csv('hpvg_and_hpvg_q5_peptide_designs_with_hcnn__pae_cutoff=5.1__seed=43__filtered_by_AF3.csv')

    ## figure out which groups the filtered-out peptides belonged to

    excluded_peptides = set(df_chosen_first['sequence'].values) - set(df_filtered['sequence'].values)

    print(excluded_peptides)

    # print out sampling_method and wt_for_design for the peptides that were excluded, i.e. from df_chosen_first
    pairs = []
    for seq in excluded_peptides:
        curr_df = df_chosen_first.loc[df_chosen_first['sequence'] == seq][['sampling_method', 'wt_for_design']]
        pairs.append((curr_df['sampling_method'].values[0], curr_df['wt_for_design'].values[0]))
    pairs = sorted(pairs)
    for pair in pairs:
        print(pair)
    
    '''
    3 from (fixed_structure, hpvg)
    3 from (fixed_structure, hpvg_q5)
    5 from (sim_anneal, hpvg)

    So, sample the following:
        - 10 from (fixed_structure, hpvg)
        - 10 from (fixed_structure, hpvg_q5)
        - 20 from (sim_anneal, hpvg)

    '''

    def sample(df, num, peptides_to_exclude):
        indices = np.random.choice(len(df), size=num, replace=False)
        sequences = list(df['peptide'].values[indices])
        while len(set(sequences).intersection(set(peptides_to_exclude))) > 0 or len(set(sequences)) < len(sequences):
            indices = np.random.choice(len(df), size=num, replace=False)
            sequences = list(df['peptide'].values[indices])
        paes = list(df['pmhc_tcr_pae'].values[indices])
        model_types = list(df['hcnn_model'].values[indices])
        return sequences, paes, model_types

    sequences_off_limits = list(df_chosen_first['sequence'].values) # initialize with peptides that were chosen first

    sequences_hpvg_fixed_structure, paes_hpvg_fixed_structure, model_types_hpvg_fixed_structure = sample(df_hpvg_fixed_structure, 10, sequences_off_limits)
    sequences_off_limits += sequences_hpvg_fixed_structure

    sequences_hpvg_q5_fixed_structure, paes_hpvg_q5_fixed_structure, model_types_hpvg_q5_fixed_structure = sample(df_hpvg_q5_fixed_structure, 10, sequences_off_limits)
    sequences_off_limits += sequences_hpvg_q5_fixed_structure

    sequences_hpvg_sim_anneal, paes_hpvg_sim_anneal, model_types_hpvg_sim_anneal = sample(df_hpvg_sim_anneal, 20, sequences_off_limits)
    sequences_off_limits += sequences_hpvg_sim_anneal

    all_sequences = sequences_hpvg_fixed_structure + sequences_hpvg_sim_anneal + sequences_hpvg_q5_fixed_structure
    all_paes = paes_hpvg_fixed_structure + paes_hpvg_sim_anneal + paes_hpvg_q5_fixed_structure
    all_model_types = model_types_hpvg_fixed_structure + model_types_hpvg_sim_anneal + model_types_hpvg_q5_fixed_structure

    assert set(all_sequences).isdisjoint(set(list(df_chosen_first['sequence'].values)))


    ## save peptides
    wt_for_design = []
    wt_for_design += ['hpvg'] * len(sequences_hpvg_fixed_structure)
    wt_for_design += ['hpvg'] * len(sequences_hpvg_sim_anneal)
    wt_for_design += ['hpvg_q5'] * len(sequences_hpvg_q5_fixed_structure)

    sampling_method = []
    sampling_method += ['fixed_structure'] * len(sequences_hpvg_fixed_structure)
    sampling_method += ['sim_anneal'] * len(sequences_hpvg_sim_anneal)
    sampling_method += ['fixed_structure'] * len(sequences_hpvg_q5_fixed_structure)

    hamming_distances_from_hpvg = [hamming_distance(seq, HPVG_SEQ) for seq in all_sequences]
    hamming_distances_from_hpvg_q5 = [hamming_distance(seq, HPVG_Q5_SEQ) for seq in all_sequences]

    df = pd.DataFrame({
        'sequence': all_sequences,
        'model': all_model_types,
        'sampling_method': sampling_method,
        'wt_for_design': wt_for_design,
        'tcrdock_pae': all_paes,
        'hamming_distances_from_hpvg': hamming_distances_from_hpvg,
        'hamming_distance_from_hpvg_q5': hamming_distances_from_hpvg_q5
    })

    df.to_csv(f'hpvg_and_hpvg_q5_peptide_designs_with_hcnn__pae_cutoff={PAE_CUTOFF}__seed={SEED}__EXTRA_FOR_AF3_TESTING.csv', index=None)

    ## make them into json files for AF3

    import json
    from copy import deepcopy

    dict_template = {
    "name": "EBV_{sequence}_1",
    "modelSeeds": [],
    "sequences": [
        {
        "proteinChain": {
            "sequence": "GSHSMRYFYTAMSRPGRGEPRFIAVGYVDDTQFVRFDSDAASPRTEPRAPWIEQEGPEYWDRNTQIFKTNTQTYRESLRNLRGYYNQSEAGSHIIQRMYGCDLGPDGRLLRGHDQSAYDGKDYIALNEDLSSWTAADTAAQITQRKWEAARVAEQLRAYLEGLCVEWLRRYLENGKETLQR",
            "count": 1
        }
        },
        {
        "proteinChain": {
            "sequence": "QVTQSPEALRLQEGESSSLNCSYTVSGLRGLFWYRQDPGKGPEFLFTLYSAGEEKEKERLKATLTKKESFLHITAPKPEDSATYLCAVQDLGTSGSRLTFGEGTQLTVNPN",
            "count": 1
        }
        },
        {
        "proteinChain": {
            "sequence": "DSGVTQTPKHLITATGQRVTLRCSPRSGDLSVYWYQQSLDQGLQFLIQYYNGEERAKGNILERFSAQQFPDLHSELNLSSLELGDSALYFCASSARSGELFFGEGSRLTVLED",
            "count": 1
        }
        }
    ]
    }

    num_sequences_per_json_file = 20

    all_peptide_sequences = all_sequences

    num_json_files = int(np.ceil(len(all_peptide_sequences) / num_sequences_per_json_file))

    for i in range(num_json_files):

        curr_json_data = []

        curr_sequences = all_peptide_sequences[i*num_sequences_per_json_file:(i+1)*num_sequences_per_json_file]

        for sequence in curr_sequences:
            curr_dict = deepcopy(dict_template)
            curr_dict['name'] = curr_dict['name'].format(sequence=sequence)
            curr_dict['sequences'].append({'proteinChain': {'sequence': sequence, 'count': 1}})
            # put the fourth sequence in second place
            curr_dict['sequences'] = [curr_dict['sequences'][0], curr_dict['sequences'][3], curr_dict['sequences'][1], curr_dict['sequences'][2]]
            curr_json_data.append(curr_dict)

        with open(f'ebv_peptide_designs_for_refinement_{i}.json', 'w') as f:
            json.dump(curr_json_data, f, indent=2)

    







