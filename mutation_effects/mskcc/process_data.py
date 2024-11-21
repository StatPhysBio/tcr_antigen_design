
import os
import numpy as np
import pandas as pd
import h5py

# I sort of got it

# mutant,is_wt,-EC50,wt_pdb,mt_pdb,mutant_chain

tcr_to_pdb_info = {
    '1': {'wt_pdb': '5d2n-filtered', 'mutant_chain': 'I', 'offset': 0},
    '2': {'wt_pdb': 'TCR2_T00000_A0201_NLVPMVATV_0_model_1_model_2_ptm_ft4', 'mutant_chain': 'A', 'offset': 375},
    '3': {'wt_pdb': 'TCR3_T00000_A0201_NLVPMVATV_0_model_1_model_2_ptm_ft4', 'mutant_chain': 'A', 'offset': 375},
    '4': {'wt_pdb': 'TCR4_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4', 'mutant_chain': 'A', 'offset': 375},
    '5': {'wt_pdb': 'TCR5_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4', 'mutant_chain': 'A', 'offset': 375},
    '6': {'wt_pdb': 'TCR6_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4', 'mutant_chain': 'A', 'offset': 375},
    '7': {'wt_pdb': 'TCR7_T00000_B2705_GRLKALCQR_0_model_1_model_2_ptm_ft4', 'mutant_chain': 'A', 'offset': 375},
}



if __name__ == '__main__':

    with h5py.File('mskcc_antagonism_fc_predictions_corrected_revised (1).h5', 'r') as f:

        data = f['EC50_fits']

        wt_peptide_idxs = data['axis0_label0'][:].astype(int)
        tcr_idxs = data['axis1_label1'][:].astype(int)
        mutant_idxs = data['axis1_label2'][:].astype(int)
    
        wt_peptide_names = [x.decode('utf-8') for x in data['axis0_level0']]
        tcr_names = [x.decode('utf-8') for x in data['axis1_level1']]
        mutant_names = [x.decode('utf-8') for x in data['axis1_level2']]

        column_1_idxs = data['block0_items_label0'][:].astype(int)
        column_2_idxs = data['block0_items_label1'][:].astype(int)
        column_1_names = [x.decode('utf-8') for x in data['block0_items_level0']]
        column_2_names = [x.decode('utf-8') for x in data['block0_items_level1']]

        chosen_column_1 = 'MAP'
        chosen_column_2 = 'log_ec50_M'
        chosen_column_1_idx = column_1_names.index(chosen_column_1)
        chosen_column_2_idx = column_2_names.index(chosen_column_2)

        mask_column_1 = column_1_idxs == chosen_column_1_idx
        mask_column_2 = column_2_idxs == chosen_column_2_idx
        mask = np.logical_and(mask_column_1, mask_column_2)
        assert np.sum(mask) == 1
        chosen_column_idx = np.arange(mask.size)[mask][0]

        for tcr_name in tcr_names:
            tcr_idx = tcr_names.index(tcr_name)
            tcr_info = tcr_to_pdb_info[tcr_name]

            curr_tcr_idxs = tcr_idxs == tcr_idx

            curr_mutant_idxs = mutant_idxs[curr_tcr_idxs]

            curr_mutants = [mutant_names[x] for x in curr_mutant_idxs]

            # convert name of WT into {wt}1{wt}
            mutants_at_pos_1 = [x for x in curr_mutants if x[1] == '1']
            wt_aa = mutants_at_pos_1[0][0]
            location_of_WT = curr_mutants.index('WT')
            curr_mutants[location_of_WT] = f'{wt_aa}1{wt_aa}'

            # offset the mutant resnums
            curr_mutant_resnums = [int(x[1:-1]) for x in curr_mutants]
            curr_mutant_resnums = [x + tcr_info['offset'] for x in curr_mutant_resnums]
            curr_mutants = [f'{x[0]}{y}{x[2]}' for x, y in zip(curr_mutants, curr_mutant_resnums)]
            wt_mutant = curr_mutants[location_of_WT]

            ec50_values = data['block0_values'][curr_tcr_idxs, chosen_column_idx]

            # shift the values by the WT value
            wt_value = ec50_values[location_of_WT]
            ec50_values = ec50_values - wt_value

            # invert the values so that higher is better
            ec50_values = -ec50_values


            # finally, make dataframe!
            df = pd.DataFrame({
                'mutant': curr_mutants,
                'is_wt': [x == wt_mutant for x in curr_mutants],
                '- delta ' + chosen_column_2: ec50_values,
                'wt_pdb': tcr_info['wt_pdb'],
                'mutant_chain': tcr_info['mutant_chain'],
            })
            df.to_csv(f'mskcc_tcr{tcr_name}_ec50_sat_mut.csv')



                          
