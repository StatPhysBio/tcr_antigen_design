
import os
import numpy as np
import pandas as pd
import h5py

tcr_to_pdb_info = {
    '8': {'wt_pdb': '6uk4_filtered', 'mutant_chain': 'C', 'offset': 0},
}


if __name__ == '__main__':

    with h5py.File('hhatv4_antagonism_fc_predictions_corrected_revised (1).h5', 'r') as f:
        # print()
        # print(f.keys())
        # print()
        # print(f['EC50_fits'].keys())
        # print()
        # print(f['EC50_fits']['axis0_label0'][:])
        # print(f['EC50_fits']['axis0_label1'][:])
        # print(f['EC50_fits']['axis0_level0'][:])
        # print(f['EC50_fits']['axis0_level1'][:])
        # print()
        # print(f['EC50_fits']['axis1_label0'][:])
        # print(f['EC50_fits']['axis1_label1'][:])
        # print(f['EC50_fits']['axis1_label2'][:])
        # print(f['EC50_fits']['axis1_level0'][:])
        # print(f['EC50_fits']['axis1_level1'][:])
        # print(f['EC50_fits']['axis1_level2'][:].shape)
        # print()
        # print(f['EC50_fits']['block0_items_label0'][:])
        # print(f['EC50_fits']['block0_items_label1'][:])
        # print(f['EC50_fits']['block0_items_level0'][:])
        # print(f['EC50_fits']['block0_items_level1'][:])
        # print(f['EC50_fits']['block0_values'][:])
        # print(f['EC50_fits']['block0_values'][:].shape)

        # exit(1)

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

            # drop the "DMSO" and "PMA/iono" mutants, and drop the corresponding idx in curr_tcr_idxs
            mutants_to_drop = ['DMSO', 'PMA/iono', 'WT', 'WildType', 'L8F']
            idxs_to_drop = []
            for mutant_to_drop in mutants_to_drop:
                if mutant_to_drop not in curr_mutants:
                    continue
                idx_to_drop = curr_mutants.index(mutant_to_drop)
                idxs_to_drop.append(idx_to_drop)
            curr_mutants = [x for x in curr_mutants if x not in mutants_to_drop]

            # convert name of WT into {wt}1{wt}
            mutants_at_pos_1 = [x for x in curr_mutants if x[1] == '1']
            wt_aa = mutants_at_pos_1[0][0]
            location_of_WT = curr_mutants.index('p8F') # the "wildtype" is the neoantigen in this case
            curr_mutants[location_of_WT] = f'{wt_aa}1{wt_aa}'

            # offset the mutant resnums
            curr_mutant_resnums = [int(x[1:-1]) for x in curr_mutants]
            curr_mutant_resnums = [x + tcr_info['offset'] for x in curr_mutant_resnums]
            curr_mutants = [f'{x[0]}{y}{x[2]}' for x, y in zip(curr_mutants, curr_mutant_resnums)]
            wt_mutant = curr_mutants[location_of_WT]

            ec50_values = data['block0_values'][curr_tcr_idxs, chosen_column_idx]
            # drop the idxs
            ec50_values = np.delete(ec50_values, idxs_to_drop)

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
            df.to_csv(f'hhatv4_tcr{tcr_name}_ec50_sat_mut.csv')


