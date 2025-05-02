
import os
import numpy as np
import pandas as pd

from tqdm import tqdm

def read_nettcr_query_file(filepath):
    # no header, but in the following order: peptide,A1,A2,A3,B1,B2,B3
    # peptides differ by rows, but the TCR info doesn't

    peptide_list = []
    with open(filepath, 'r') as f:
        for line in f:
            peptide,A1,A2,A3,B1,B2,B3 = line.strip('\n').split(',')
            peptide_list.append(peptide)
    
    return {
        'A1': A1,
        'A2': A2,
        'A3': A3,
        'B1': B1,
        'B2': B2,
        'B3': B3,
        'peptide_list': peptide_list
    }

def make_presence_bool_df(df_nettcr, test_data, columns_to_consider):

    df_presence = pd.DataFrame(columns=columns_to_consider)
    for i, row in tqdm(df_nettcr.iterrows(), total=len(df_nettcr)):
        new_row = {}
        for col in columns_to_consider:
            if pd.isna(row[col]):
                new_row[col] = np.nan
            else:
                if col == 'peptide':
                    new_row[col] = row[col] in test_data['peptide_list']
                else:
                    new_row[col] = row[col] == test_data[col]
        
        df_presence = pd.concat([df_presence, pd.DataFrame(new_row, index=[0])], ignore_index=True)

    return df_presence


def print_presence_and_matches(df_nettcr, df_presence):

    print(np.nansum(df_presence['A1'].values))
    print(np.nansum(df_presence['A2'].values))
    print(np.nansum(df_presence['A3'].values))
    print(np.nansum(df_presence['B1'].values))
    print(np.nansum(df_presence['B2'].values))
    print(np.nansum(df_presence['B3'].values))
    print(np.nansum(df_presence['peptide'].values))

    num_matches_per_row = np.nansum(df_presence.values, axis=1)
    print(np.where(num_matches_per_row == 3)[0])
    print(np.where(num_matches_per_row == 4)[0])
    print(np.where(num_matches_per_row == 5)[0])
    print(np.where(num_matches_per_row == 6)[0])
    print(np.where(num_matches_per_row == 7)[0])


df_nettcr = pd.read_csv('../../../../Downloads/nettcr_2_2_full_dataset.csv')

columns_to_consider = 'A1,A2,A3,B1,B2,B3,peptide'.split(',')


# print('nyeso')
# test_data = read_nettcr_query_file('nyeso/nettcr_input_nyeso.txt')
# df_presence = make_presence_bool_df(df_nettcr, test_data, columns_to_consider)
# print_presence_and_matches(df_nettcr, df_presence)
# print()
# print()

# print('tax')
# test_data = read_nettcr_query_file('tax/nettcr_input_tax.txt')
# df_presence = make_presence_bool_df(df_nettcr, test_data, columns_to_consider)
# print_presence_and_matches(df_nettcr, df_presence)
# print()
# print()

# print('tcr1')
# test_data = read_nettcr_query_file('mskcc/nettcr_input_tcr1.txt')
# df_presence = make_presence_bool_df(df_nettcr, test_data, columns_to_consider)
# print_presence_and_matches(df_nettcr, df_presence)
# print()
# print()

# print('tcr2')
# test_data = read_nettcr_query_file('mskcc/nettcr_input_tcr2.txt')
# df_presence = make_presence_bool_df(df_nettcr, test_data, columns_to_consider)
# print_presence_and_matches(df_nettcr, df_presence)
# print()
# print()

# print('tcr3')
# test_data = read_nettcr_query_file('mskcc/nettcr_input_tcr3.txt')
# df_presence = make_presence_bool_df(df_nettcr, test_data, columns_to_consider)
# print_presence_and_matches(df_nettcr, df_presence)
# print()
# print()

# print('tcr4')
# test_data = read_nettcr_query_file('mskcc/nettcr_input_tcr4.txt')
# df_presence = make_presence_bool_df(df_nettcr, test_data, columns_to_consider)
# print_presence_and_matches(df_nettcr, df_presence)
# print()
# print()

# print('tcr5')
# test_data = read_nettcr_query_file('mskcc/nettcr_input_tcr5.txt')
# df_presence = make_presence_bool_df(df_nettcr, test_data, columns_to_consider)
# print_presence_and_matches(df_nettcr, df_presence)
# print()
# print()

# print('tcr6')
# test_data = read_nettcr_query_file('mskcc/nettcr_input_tcr6.txt')
# df_presence = make_presence_bool_df(df_nettcr, test_data, columns_to_consider)
# print_presence_and_matches(df_nettcr, df_presence)
# print()
# print()

print('tcr7')
test_data = read_nettcr_query_file('mskcc/nettcr_input_tcr7.txt')
df_presence = make_presence_bool_df(df_nettcr, test_data, columns_to_consider)
print_presence_and_matches(df_nettcr, df_presence)
print()
print()

