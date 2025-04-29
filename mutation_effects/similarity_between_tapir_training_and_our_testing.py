
import os
import numpy as np
import pandas as pd

from tqdm import tqdm

df_tapir = pd.read_excel('../../../../Downloads/media-2.xlsx', sheet_name='Supp Table1')

# remove header, make first row header
df_tapir.columns = df_tapir.iloc[0]
df_tapir = df_tapir[1:]

# rename columns
df_tapir = df_tapir.rename(columns={'alpha v': 'alpha_v',
                                  'alpha j': 'alpha_j',
                                  'alpha cdr3': 'alpha_cdr3',
                                  'beta v': 'beta_v',
                                  'beta j': 'beta_j',
                                  'beta cdr3': 'beta_cdr3'})


# df_nyeso = pd.read_csv('nyeso/nyeso_tapir.tsv', sep='\t')
# mhc_nyeso = 'HLA-A*02'
# peptides_nyeso = pd.read_csv('nyeso/nyeso_peptide_kd_closest.csv')['sequence'].values

# # make a new df like df_tapir, where the columns are the same as df_tapir, and the values are whether that is a match with df_nyeso, except if the df_tapir value is NaN
# df_nyeso_tapir = pd.DataFrame(columns=df_tapir.columns)
# for i, row in tqdm(df_tapir.iterrows(), total=df_tapir.shape[0]):
#     new_row = {}
#     for col in df_tapir.columns:
#         if pd.isna(row[col]):
#             new_row[col] = np.nan
#         else:
#             if col == 'antigen':
#                 new_row[col] = row[col] in peptides_nyeso
#             elif col == 'mhc':
#                 new_row[col] = row[col] == mhc_nyeso
#             else:
#                 new_row[col] = row[col] in df_nyeso[col].values
    
#     df_nyeso_tapir = pd.concat([df_nyeso_tapir, pd.DataFrame(new_row, index=[0])], ignore_index=True)

# print(df_nyeso_tapir)

# df_nyeso_tapir.to_csv('nyeso_tapir_matches.csv', index=False)

# df_nyeso_tapir = pd.read_csv('nyeso_tapir_matches.csv')

# print(np.nansum(df_nyeso_tapir['alpha_v'].values))
# print(np.nansum(df_nyeso_tapir['alpha_j'].values))
# print(np.nansum(df_nyeso_tapir['alpha_cdr3'].values))
# print(np.nansum(df_nyeso_tapir['beta_v'].values))
# print(np.nansum(df_nyeso_tapir['beta_j'].values))
# print(np.nansum(df_nyeso_tapir['beta_cdr3'].values))
# print(np.nansum(df_nyeso_tapir['antigen'].values))
# print(np.nansum(df_nyeso_tapir['mhc'].values))

# num_matches_per_row = np.nansum(df_nyeso_tapir.values, axis=1)
# print(np.where(num_matches_per_row == 3)[0])
# print(np.where(num_matches_per_row == 4)[0])
# print(np.where(num_matches_per_row == 5)[0])
# print(np.where(num_matches_per_row == 6)[0])
# print(np.where(num_matches_per_row == 7)[0])
# print(np.where(num_matches_per_row == 8)[0])


# df_tax = pd.read_csv('tax/tax_tapir.tsv', sep='\t')
# mhc_tax = 'HLA-A*02'
# peptides_tax = pd.read_csv('tax/tax_peptide_kd_closest.csv')['sequence'].values

# # df_tax_tapir = pd.DataFrame(columns=df_tapir.columns)
# # for i, row in tqdm(df_tapir.iterrows(), total=df_tapir.shape[0]):
# #     new_row = {}
# #     for col in df_tapir.columns:
# #         if pd.isna(row[col]):
# #             new_row[col] = np.nan
# #         else:
# #             if col == 'antigen':
# #                 new_row[col] = row[col] in peptides_tax
# #             elif col == 'mhc':
# #                 new_row[col] = row[col] == mhc_tax
# #             else:
# #                 new_row[col] = row[col] in df_tax[col].values
    
# #     df_tax_tapir = pd.concat([df_tax_tapir, pd.DataFrame(new_row, index=[0])], ignore_index=True)

# # print(df_tax_tapir)

# # df_tax_tapir.to_csv('tax_tapir_matches.csv', index=False)

# df_tax_tapir = pd.read_csv('tax_tapir_matches.csv')

# print(np.nansum(df_tax_tapir['alpha_v'].values))
# print(np.nansum(df_tax_tapir['alpha_j'].values))
# print(np.nansum(df_tax_tapir['alpha_cdr3'].values))
# print(np.nansum(df_tax_tapir['beta_v'].values))
# print(np.nansum(df_tax_tapir['beta_j'].values))
# print(np.nansum(df_tax_tapir['beta_cdr3'].values))
# print(np.nansum(df_tax_tapir['antigen'].values))
# print(np.nansum(df_tax_tapir['mhc'].values))

# num_matches_per_row = np.nansum(df_tax_tapir.values, axis=1)
# print(np.where(num_matches_per_row == 3)[0])
# print(np.where(num_matches_per_row == 4)[0])
# print(np.where(num_matches_per_row == 5)[0])
# print(np.where(num_matches_per_row == 6)[0])
# print(np.where(num_matches_per_row == 7)[0])
# print(np.where(num_matches_per_row == 8)[0])

# print(df_tax_tapir.loc[np.where(num_matches_per_row == 7)[0]])
# print(df_tapir.loc[np.where(num_matches_per_row == 7)[0]])


# df_tcr1 = pd.read_csv('mskcc/mskcc_tcr1_tapir.tsv', sep='\t')
# mhc_tcr1 = 'HLA-A*02'
# peptides_tcr1 = pd.read_csv('mskcc/mskcc_tcr1_ec50_sat_mut_af3.csv')['sequence'].values

# # df_tcr1_tapir = pd.DataFrame(columns=df_tapir.columns)
# # for i, row in tqdm(df_tapir.iterrows(), total=df_tapir.shape[0]):
# #     new_row = {}
# #     for col in df_tapir.columns:
# #         if pd.isna(row[col]):
# #             new_row[col] = np.nan
# #         else:
# #             if col == 'antigen':
# #                 new_row[col] = row[col] in peptides_tcr1
# #             elif col == 'mhc':
# #                 new_row[col] = row[col] == mhc_tcr1
# #             else:
# #                 new_row[col] = row[col] in df_tcr1[col].values
    
# #     df_tcr1_tapir = pd.concat([df_tcr1_tapir, pd.DataFrame(new_row, index=[0])], ignore_index=True)

# # print(df_tcr1_tapir)

# # df_tcr1_tapir.to_csv('tcr1_tapir_matches.csv', index=False)

# df_tcr1_tapir = pd.read_csv('tcr1_tapir_matches.csv')

# print(np.nansum(df_tcr1_tapir['alpha_v'].values))
# print(np.nansum(df_tcr1_tapir['alpha_j'].values))
# print(np.nansum(df_tcr1_tapir['alpha_cdr3'].values))
# print(np.nansum(df_tcr1_tapir['beta_v'].values))
# print(np.nansum(df_tcr1_tapir['beta_j'].values))
# print(np.nansum(df_tcr1_tapir['beta_cdr3'].values))
# print(np.nansum(df_tcr1_tapir['antigen'].values))
# print(np.nansum(df_tcr1_tapir['mhc'].values))

# num_matches_per_row = np.nansum(df_tcr1_tapir.values, axis=1)
# print(np.where(num_matches_per_row == 3)[0])
# print(np.where(num_matches_per_row == 4)[0])
# print(np.where(num_matches_per_row == 5)[0])
# print(np.where(num_matches_per_row == 6)[0])
# print(np.where(num_matches_per_row == 7)[0])
# print(np.where(num_matches_per_row == 8)[0])

# print(df_tapir.iloc[np.where(num_matches_per_row == 8)[0]])



# df_tcr2 = pd.read_csv('mskcc/mskcc_tcr2_tapir.tsv', sep='\t')
# mhc_tcr2 = 'HLA-A*02'
# peptides_tcr2 = pd.read_csv('mskcc/mskcc_tcr2_ec50_sat_mut_af3.csv')['sequence'].values

# # df_tcr2_tapir = pd.DataFrame(columns=df_tapir.columns)
# # for i, row in tqdm(df_tapir.iterrows(), total=df_tapir.shape[0]):
# #     new_row = {}
# #     for col in df_tapir.columns:
# #         if pd.isna(row[col]):
# #             new_row[col] = np.nan
# #         else:
# #             if col == 'antigen':
# #                 new_row[col] = row[col] in peptides_tcr2
# #             elif col == 'mhc':
# #                 new_row[col] = row[col] == mhc_tcr2
# #             else:
# #                 new_row[col] = row[col] in df_tcr2[col].values
    
# #     df_tcr2_tapir = pd.concat([df_tcr2_tapir, pd.DataFrame(new_row, index=[0])], ignore_index=True)

# # print(df_tcr2_tapir)

# # df_tcr2_tapir.to_csv('tcr2_tapir_matches.csv', index=False)

# df_tcr2_tapir = pd.read_csv('tcr2_tapir_matches.csv')

# print(np.nansum(df_tcr2_tapir['alpha_v'].values))
# print(np.nansum(df_tcr2_tapir['alpha_j'].values))
# print(np.nansum(df_tcr2_tapir['alpha_cdr3'].values))
# print(np.nansum(df_tcr2_tapir['beta_v'].values))
# print(np.nansum(df_tcr2_tapir['beta_j'].values))
# print(np.nansum(df_tcr2_tapir['beta_cdr3'].values))
# print(np.nansum(df_tcr2_tapir['antigen'].values))
# print(np.nansum(df_tcr2_tapir['mhc'].values))

# num_matches_per_row = np.nansum(df_tcr2_tapir.values, axis=1)
# print(np.where(num_matches_per_row == 3)[0])
# print(np.where(num_matches_per_row == 4)[0])
# print(np.where(num_matches_per_row == 5)[0])
# print(np.where(num_matches_per_row == 6)[0])
# print(np.where(num_matches_per_row == 7)[0])
# print(np.where(num_matches_per_row == 8)[0])

# print(df_tcr2_tapir.iloc[np.where(num_matches_per_row == 7)[0]])
# print(df_tapir.iloc[np.where(num_matches_per_row == 7)[0]])



# df_tcr3 = pd.read_csv('mskcc/mskcc_tcr3_tapir.tsv', sep='\t')
# mhc_tcr3 = 'HLA-A*02'
# peptides_tcr3 = pd.read_csv('mskcc/mskcc_tcr3_ec50_sat_mut_af3.csv')['sequence'].values

# # df_tcr3_tapir = pd.DataFrame(columns=df_tapir.columns)
# # for i, row in tqdm(df_tapir.iterrows(), total=df_tapir.shape[0]):
# #     new_row = {}
# #     for col in df_tapir.columns:
# #         if pd.isna(row[col]):
# #             new_row[col] = np.nan
# #         else:
# #             if col == 'antigen':
# #                 new_row[col] = row[col] in peptides_tcr3
# #             elif col == 'mhc':
# #                 new_row[col] = row[col] == mhc_tcr3
# #             else:
# #                 new_row[col] = row[col] in df_tcr3[col].values
    
# #     df_tcr3_tapir = pd.concat([df_tcr3_tapir, pd.DataFrame(new_row, index=[0])], ignore_index=True)

# # print(df_tcr3_tapir)

# # df_tcr3_tapir.to_csv('tcr3_tapir_matches.csv', index=False)

# df_tcr3_tapir = pd.read_csv('tcr3_tapir_matches.csv')

# print(np.nansum(df_tcr3_tapir['alpha_v'].values))
# print(np.nansum(df_tcr3_tapir['alpha_j'].values))
# print(np.nansum(df_tcr3_tapir['alpha_cdr3'].values))
# print(np.nansum(df_tcr3_tapir['beta_v'].values))
# print(np.nansum(df_tcr3_tapir['beta_j'].values))
# print(np.nansum(df_tcr3_tapir['beta_cdr3'].values))
# print(np.nansum(df_tcr3_tapir['antigen'].values))
# print(np.nansum(df_tcr3_tapir['mhc'].values))

# num_matches_per_row = np.nansum(df_tcr3_tapir.values, axis=1)
# print(np.where(num_matches_per_row == 3)[0])
# print(np.where(num_matches_per_row == 4)[0])
# print(np.where(num_matches_per_row == 5)[0])
# print(np.where(num_matches_per_row == 6)[0])
# print(np.where(num_matches_per_row == 7)[0])
# print(np.where(num_matches_per_row == 8)[0])

# print(df_tcr3_tapir.iloc[np.where(num_matches_per_row == 5)[0]])
# print(df_tapir.iloc[np.where(num_matches_per_row == 5)[0]])



# df_tcr4 = pd.read_csv('mskcc/mskcc_tcr4_tapir.tsv', sep='\t')
# mhc_tcr4 = 'HLA-A*02'
# peptides_tcr4 = pd.read_csv('mskcc/mskcc_tcr4_ec50_sat_mut_af3.csv')['sequence'].values

# df_tcr4_tapir = pd.DataFrame(columns=df_tapir.columns)
# for i, row in tqdm(df_tapir.iterrows(), total=df_tapir.shape[0]):
#     new_row = {}
#     for col in df_tapir.columns:
#         if pd.isna(row[col]):
#             new_row[col] = np.nan
#         else:
#             if col == 'antigen':
#                 new_row[col] = row[col] in peptides_tcr4
#             elif col == 'mhc':
#                 new_row[col] = row[col] == mhc_tcr4
#             else:
#                 new_row[col] = row[col] in df_tcr4[col].values
    
#     df_tcr4_tapir = pd.concat([df_tcr4_tapir, pd.DataFrame(new_row, index=[0])], ignore_index=True)

# print(df_tcr4_tapir)

# df_tcr4_tapir.to_csv('tcr4_tapir_matches.csv', index=False)

# df_tcr4_tapir = pd.read_csv('tcr4_tapir_matches.csv')

# print(np.nansum(df_tcr4_tapir['alpha_v'].values))
# print(np.nansum(df_tcr4_tapir['alpha_j'].values))
# print(np.nansum(df_tcr4_tapir['alpha_cdr3'].values))
# print(np.nansum(df_tcr4_tapir['beta_v'].values))
# print(np.nansum(df_tcr4_tapir['beta_j'].values))
# print(np.nansum(df_tcr4_tapir['beta_cdr3'].values))
# print(np.nansum(df_tcr4_tapir['antigen'].values))
# print(np.nansum(df_tcr4_tapir['mhc'].values))

# num_matches_per_row = np.nansum(df_tcr4_tapir.values, axis=1)
# print(np.where(num_matches_per_row == 3)[0])
# print(np.where(num_matches_per_row == 4)[0])
# print(np.where(num_matches_per_row == 5)[0])
# print(np.where(num_matches_per_row == 6)[0])
# print(np.where(num_matches_per_row == 7)[0])
# print(np.where(num_matches_per_row == 8)[0])

# print(df_tcr4_tapir.iloc[np.where(num_matches_per_row == 5)[0]])
# print(df_tapir.iloc[np.where(num_matches_per_row == 5)[0]])



# df_tcr5 = pd.read_csv('mskcc/mskcc_tcr5_tapir.tsv', sep='\t')
# mhc_tcr5 = 'HLA-A*02'
# peptides_tcr5 = pd.read_csv('mskcc/mskcc_tcr5_ec50_sat_mut_af3.csv')['sequence'].values

# df_tcr5_tapir = pd.DataFrame(columns=df_tapir.columns)
# for i, row in tqdm(df_tapir.iterrows(), total=df_tapir.shape[0]):
#     new_row = {}
#     for col in df_tapir.columns:
#         if pd.isna(row[col]):
#             new_row[col] = np.nan
#         else:
#             if col == 'antigen':
#                 new_row[col] = row[col] in peptides_tcr5
#             elif col == 'mhc':
#                 new_row[col] = row[col] == mhc_tcr5
#             else:
#                 new_row[col] = row[col] in df_tcr5[col].values
    
#     df_tcr5_tapir = pd.concat([df_tcr5_tapir, pd.DataFrame(new_row, index=[0])], ignore_index=True)

# print(df_tcr5_tapir)

# df_tcr5_tapir.to_csv('tcr5_tapir_matches.csv', index=False)

# df_tcr5_tapir = pd.read_csv('tcr5_tapir_matches.csv')

# print(np.nansum(df_tcr5_tapir['alpha_v'].values))
# print(np.nansum(df_tcr5_tapir['alpha_j'].values))
# print(np.nansum(df_tcr5_tapir['alpha_cdr3'].values))
# print(np.nansum(df_tcr5_tapir['beta_v'].values))
# print(np.nansum(df_tcr5_tapir['beta_j'].values))
# print(np.nansum(df_tcr5_tapir['beta_cdr3'].values))
# print(np.nansum(df_tcr5_tapir['antigen'].values))
# print(np.nansum(df_tcr5_tapir['mhc'].values))

# num_matches_per_row = np.nansum(df_tcr5_tapir.values, axis=1)
# print(np.where(num_matches_per_row == 3)[0])
# print(np.where(num_matches_per_row == 4)[0])
# print(np.where(num_matches_per_row == 5)[0])
# print(np.where(num_matches_per_row == 6)[0])
# print(np.where(num_matches_per_row == 7)[0])
# print(np.where(num_matches_per_row == 8)[0])

# print(df_tcr5_tapir.iloc[np.where(num_matches_per_row == 8)[0]])
# print(df_tapir.iloc[np.where(num_matches_per_row == 8)[0]])


# df_tcr6 = pd.read_csv('mskcc/mskcc_tcr6_tapir.tsv', sep='\t')
# mhc_tcr6 = 'HLA-A*02'
# peptides_tcr6 = pd.read_csv('mskcc/mskcc_tcr6_ec50_sat_mut_af3.csv')['sequence'].values

# df_tcr6_tapir = pd.DataFrame(columns=df_tapir.columns)
# for i, row in tqdm(df_tapir.iterrows(), total=df_tapir.shape[0]):
#     new_row = {}
#     for col in df_tapir.columns:
#         if pd.isna(row[col]):
#             new_row[col] = np.nan
#         else:
#             if col == 'antigen':
#                 new_row[col] = row[col] in peptides_tcr6
#             elif col == 'mhc':
#                 new_row[col] = row[col] == mhc_tcr6
#             else:
#                 new_row[col] = row[col] in df_tcr6[col].values
    
#     df_tcr6_tapir = pd.concat([df_tcr6_tapir, pd.DataFrame(new_row, index=[0])], ignore_index=True)

# print(df_tcr6_tapir)

# df_tcr6_tapir.to_csv('tcr6_tapir_matches.csv', index=False)

# df_tcr6_tapir = pd.read_csv('tcr6_tapir_matches.csv')

# print(np.nansum(df_tcr6_tapir['alpha_v'].values))
# print(np.nansum(df_tcr6_tapir['alpha_j'].values))
# print(np.nansum(df_tcr6_tapir['alpha_cdr3'].values))
# print(np.nansum(df_tcr6_tapir['beta_v'].values))
# print(np.nansum(df_tcr6_tapir['beta_j'].values))
# print(np.nansum(df_tcr6_tapir['beta_cdr3'].values))
# print(np.nansum(df_tcr6_tapir['antigen'].values))
# print(np.nansum(df_tcr6_tapir['mhc'].values))

# num_matches_per_row = np.nansum(df_tcr6_tapir.values, axis=1)
# print(np.where(num_matches_per_row == 3)[0])
# print(np.where(num_matches_per_row == 4)[0])
# print(np.where(num_matches_per_row == 5)[0])
# print(np.where(num_matches_per_row == 6)[0])
# print(np.where(num_matches_per_row == 7)[0])
# print(np.where(num_matches_per_row == 8)[0])

# print(df_tcr6_tapir.iloc[np.where(num_matches_per_row == 8)[0]])
# print(df_tapir.iloc[np.where(num_matches_per_row == 8)[0]])



# df_tcr7 = pd.read_csv('mskcc/mskcc_tcr7_tapir.tsv', sep='\t')
# mhc_tcr7 = 'HLA-B*27'
# peptides_tcr7 = pd.read_csv('mskcc/mskcc_tcr7_ec50_sat_mut_af3.csv')['sequence'].values

# df_tcr7_tapir = pd.DataFrame(columns=df_tapir.columns)
# for i, row in tqdm(df_tapir.iterrows(), total=df_tapir.shape[0]):
#     new_row = {}
#     for col in df_tapir.columns:
#         if pd.isna(row[col]):
#             new_row[col] = np.nan
#         else:
#             if col == 'antigen':
#                 new_row[col] = row[col] in peptides_tcr7
#             elif col == 'mhc':
#                 new_row[col] = row[col] == mhc_tcr7
#             else:
#                 new_row[col] = row[col] in df_tcr7[col].values
    
#     df_tcr7_tapir = pd.concat([df_tcr7_tapir, pd.DataFrame(new_row, index=[0])], ignore_index=True)

# print(df_tcr7_tapir)

# df_tcr7_tapir.to_csv('tcr7_tapir_matches.csv', index=False)

df_tcr7_tapir = pd.read_csv('tcr7_tapir_matches.csv')

print(np.nansum(df_tcr7_tapir['alpha_v'].values))
print(np.nansum(df_tcr7_tapir['alpha_j'].values))
print(np.nansum(df_tcr7_tapir['alpha_cdr3'].values))
print(np.nansum(df_tcr7_tapir['beta_v'].values))
print(np.nansum(df_tcr7_tapir['beta_j'].values))
print(np.nansum(df_tcr7_tapir['beta_cdr3'].values))
print(np.nansum(df_tcr7_tapir['antigen'].values))
print(np.nansum(df_tcr7_tapir['mhc'].values))

num_matches_per_row = np.nansum(df_tcr7_tapir.values, axis=1)
print(np.where(num_matches_per_row == 3)[0])
print(np.where(num_matches_per_row == 4)[0])
print(np.where(num_matches_per_row == 5)[0])
print(np.where(num_matches_per_row == 6)[0])
print(np.where(num_matches_per_row == 7)[0])
print(np.where(num_matches_per_row == 8)[0])

print(df_tcr7_tapir.iloc[np.where(num_matches_per_row == 3)[0]])
print(df_tapir.iloc[np.where(num_matches_per_row == 3)[0]])

