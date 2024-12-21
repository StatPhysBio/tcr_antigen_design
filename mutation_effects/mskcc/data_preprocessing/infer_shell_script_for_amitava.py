#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 21:32:05 2023

@author: zacharysethna
"""

import os
import numpy as np
import time
import sys
# sys.path.insert(0, '/Users/zacharysethna/MSKresearch/crossreactivity/')
from InferEpitopeDist import TCRpMHCAvidity, InferEpitopeDist

#%
import json
def save_inferred_model(av_m, m, outfile_name):
    outmodel_dict = {}
    for prop, val in m.items():
        if type(val) is float:
            outmodel_dict[prop] = val
        elif prop == 'euclid_coords':
            outmodel_dict[prop] = {aa: list(m['euclid_coords'][i]) for i, aa in enumerate(av_m.amino_acids)}
        elif prop == 'M_ab':
            outmodel_dict[prop] = {aaA + '->' + aaB: m['M_ab'][i, j] for i, aaA in enumerate(av_m.amino_acids) for j, aaB in enumerate(av_m.amino_acids)}
        elif type(val) is np.ndarray:
            outmodel_dict[prop] = list(val)
                    
#    outmodel_dict = {#'blosum62_reg': m['blosum62_reg'],
#                     'log_units': m['log_units'],
#                     'd_i': list(m['d_i']),
#                     #'euclid_coords': {aa: list(m['euclid_coords'][i]) for i, aa in enumerate(av_m.amino_acids)},
#                     'M_ab': {aaA + '->' + aaB: m['M_ab'][i, j] for i, aaA in enumerate(av_m.amino_acids) for j, aaB in enumerate(av_m.amino_acids)}}
    with open(outfile_name, 'w') as outfile:
        json.dump(outmodel_dict, outfile)

data_folder = './'

cmv_peptide = 'NLVPMVATV'
gp100_peptide = 'IMDQVPFSV'
neopeptide = 'GRLKALCQR'
tcr_pMHC_pairs = {
    'NLV_tcr1': cmv_peptide,
    'NLV_tcr2': cmv_peptide,
    'NLV_tcr3': cmv_peptide,
    'IMD_tcr4': gp100_peptide,
    'IMD_tcr5': gp100_peptide,
    'IMD_tcr6': gp100_peptide,
    'GRL_tcr7': neopeptide
}


start_time = time.time()
K_a_reg = 0.001
log_Ka_range = [-6, 6] # original: [-4, 4]
all_tcr_info = {}


# make the reactivity curves tsv files on the fly, extracting them from avidity_curves.xlsx
import pandas as pd

avidity_curves = pd.read_excel(os.path.join(data_folder, 'avidity_curves.xlsx')).astype(str)
columns_dict = {
    'NLV_tcr1': [0, 1, 2, 3],
    'NLV_tcr2': [5, 6, 7, 8],
    'NLV_tcr3': [10, 11, 12, 13],
    'IMD_tcr4': [15, 16, 17, 18],
    'IMD_tcr5': [20, 21, 22, 23],
    'IMD_tcr6': [25, 26, 27, 28],
    'GRL_tcr7': [30, 31, 32, 33]
}

for tcr, columns in columns_dict.items():
    tcr_avidity = avidity_curves.iloc[2:, columns]
    tcr_avidity.columns = avidity_curves.iloc[1, columns]
    tcr_avidity.to_csv(os.path.join(data_folder, tcr + '_reactivity_curves.tsv'), sep='\t', index=False)


# no wt datafiles loaded
for tcr in tcr_pMHC_pairs.keys():
    all_tcr_info[tcr] = TCRpMHCAvidity(mut_pep_avidity_infile = os.path.join(data_folder, tcr + '_reactivity_curves.tsv'), K_a_reg = K_a_reg, log_Ka_range= log_Ka_range)

print('Loaded all tcrs in %.2f seconds'%(time.time() - start_time))

for tcr in tcr_pMHC_pairs:
    ec50s = all_tcr_info[tcr].return_all_obs_C()
    with open(tcr + '_ec50s.json', 'w+') as f:
        json.dump(ec50s, f, indent=4)


# infer_ep_dist = InferEpitopeDist()
# #%
# nlv_tcrs = {tcr_name: c_tcr for tcr_name, c_tcr in all_tcr_info.items() if tcr_name.startswith('NLV')}

# #%%
# infer_ep_dist.n_dim = 2
# infer_ep_dist.reactivity_cutoff = 0.1
# infer_ep_dist.blosum62_reg = 0

# print('Starting NLV model')
# nlv_model_w_reg = infer_ep_dist.infer_euclid_coord_model_from_tcrs(list(nlv_tcrs.values()), seed = 17, use_all_pep_combos = True)
# #%%
# outdir = './'
# save_inferred_model(infer_ep_dist, nlv_model_w_reg, os.path.join(outdir, 'nlv_model_w_reg_all_combos_model_bl62reg_%.0e_rt_%.0e.json'%(infer_ep_dist.blosum62_reg,infer_ep_dist.reactivity_cutoff)))


