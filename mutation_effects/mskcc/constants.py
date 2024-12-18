
# PDB_TO_PEP_INFO = {
#     '5d2n-filtered': ['I', 1, 'NLVPMVATV'],
#     'TCR2_T00000_A0201_NLVPMVATV_0_model_1_model_2_ptm_ft4': ['A', 376, 'NLVPMVATV'],
#     'TCR3_T00000_A0201_NLVPMVATV_0_model_1_model_2_ptm_ft4': ['A', 376, 'NLVPMVATV'],
#     'TCR4_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4': ['A', 376, 'IMDQVPFSV'],
#     'TCR5_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4': ['A', 376, 'IMDQVPFSV'],
#     'TCR6_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4': ['A', 376, 'IMDQVPFSV'],
#     'TCR7_T00000_B2705_GRLKALCQR_0_model_1_model_2_ptm_ft4': ['A', 376, 'GRLKALCQR']
# }

# TCR_TO_PDB = {
#     '1': '5d2n-filtered',
#     '2': 'TCR2_T00000_A0201_NLVPMVATV_0_model_1_model_2_ptm_ft4',
#     '3': 'TCR3_T00000_A0201_NLVPMVATV_0_model_1_model_2_ptm_ft4',
#     '4': 'TCR4_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4',
#     '5': 'TCR5_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4',
#     '6': 'TCR6_T00000_A0201_IMDQVPFSV_0_model_1_model_2_ptm_ft4',
#     '7': 'TCR7_T00000_B2705_GRLKALCQR_0_model_1_model_2_ptm_ft4'
# }

PDB_TO_PEP_INFO = {
    '5d2n-filtered': ['I', 1, 'NLVPMVATV'],
    'af3_tcr2': ['B', 1, 'NLVPMVATV'],
    'af3_tcr3': ['B', 1, 'NLVPMVATV'],
    'af3_tcr4': ['B', 1, 'IMDQVPFSV'],
    'af3_tcr5': ['B', 1, 'IMDQVPFSV'],
    'af3_tcr6': ['B', 1, 'IMDQVPFSV'],
    'af3_tcr7': ['B', 1, 'GRLKALCQR']
}

TCR_TO_PDB = {
    '1': '5d2n-filtered',
    '2': 'af3_tcr2',
    '3': 'af3_tcr3',
    '4': 'af3_tcr4',
    '5': 'af3_tcr5',
    '6': 'af3_tcr6',
    '7': 'af3_tcr7'
}

PDBS = list(PDB_TO_PEP_INFO.keys())
