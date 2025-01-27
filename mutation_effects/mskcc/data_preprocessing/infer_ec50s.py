#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 21:32:05 2023

@author: zacharysethna

Edited Sat Dec 21 2024

@author: Armita Nourmohammad

Edited Sun Jan 26 2025

@author: Gian Marco Visani
"""

import os
import json
import numpy as np
import scipy.optimize as opt
import os
import pandas as pd


class TCRpMHCAvidity(object):
    def __init__(self, **kwarg):
        if 'wt_peptide' in kwarg:
            self.wt_peptide = kwarg['wt_peptide']
        else:
            self.wt_peptide = None
        
        if 'n_reg' in kwarg:
            self.n_reg = kwarg['n_reg']
        else:
            self.n_reg = 0.01
            
        if 'A_reg' in kwarg:
            self.A_reg = kwarg['A_reg']
        else:
            self.A_reg = 0.01
        
        if 'bgrn_reg' in kwarg:
            self.bgrn_reg = kwarg['bgrn_reg']
        else:
            self.bgrn_reg = 0.1
        
        
        if 'K_a_reg' in kwarg:
            self.K_a_reg = kwarg['K_a_reg']
        else:
            self.K_a_reg = 0
            
        if 'log_Ka_range' in kwarg:
            self.log_Ka_range = kwarg['log_Ka_range']
        else:
            self.log_Ka_range = [-4, 4]
            
        self.peptide_reactivity_data = {}
        self.peptide_avidity_fits = {}
        
        if 'mut_pep_avidity_infile' in kwarg:
            self.load_mut_pep_avidity_data(kwarg['mut_pep_avidity_infile'])
        
        self.wt_avidity = None
        
        if 'wt_avidity_infile' in kwarg and os.path.exists(kwarg['wt_avidity_infile']):
            self.load_wt_avidity_data(kwarg['wt_avidity_infile'])
        elif any([pep[0] == pep[-1] for pep in self.peptide_reactivity_data.keys()]):
            wt_pep_rep = [pep for pep in self.peptide_reactivity_data.keys() if pep[0] == pep[-1]][0]
            self.wt_reactivity = self.peptide_reactivity_data[wt_pep_rep]
            self.wt_fit = self.peptide_avidity_fits[wt_pep_rep]
        
        if 'wt_avidity' in kwarg:
            self.wt_avidity = kwarg['wt_avidity']
        else:
            try:
                self.wt_avidity = self.wt_fit['K_a']
            except:
                pass
        
    def load_wt_avidity_data(self, infile_name):
        with open(infile_name, 'r') as infile:
            all_L = [l.split('\t') for l in infile.read().strip().split('\n')]
        self.wt_reactivity = {float(l[0]): float(l[1]) for l in all_L[1:]}
        self.wt_fit = self.fit_indiv_t_react_curve(self.wt_reactivity)
    
    def load_mut_pep_avidity_data(self, infile_name):
        with open(infile_name, 'r') as infile:
            all_L = [l.split('\t') for l in infile.read().strip().split('\n')]
            
        concs = [float(c) for c in all_L[0][1:]]
        
        for l in all_L[1:]:
            c_react_dict = {}
            for c, react in zip(concs, l[1:]):
                if react in {'N/A', 'nan'}: continue
                c_react_dict[c] = float(react)
        
            if l[0] not in self.peptide_reactivity_data:
                self.peptide_reactivity_data[l[0]] = c_react_dict
            else:
                for c, react in c_react_dict.items():
                    self.peptide_reactivity_data[l[0]][c] = react
        if self.wt_peptide is None: self.determine_wt_sequence()
        
        self.fit_mut_pep_avidity()
    
    def determine_wt_sequence(self):            
        c_peptides = self.peptide_reactivity_data.keys()
        mut_pos_inds = np.array(sorted(set([int(pep[1:-1]) for pep in c_peptides])))
        try:
            if all(mut_pos_inds == np.array(range(1, len(mut_pos_inds) +1))):
                wt_aa_by_pos = [set() for _ in mut_pos_inds]
                for pep in c_peptides:
                    wt_aa_by_pos[int(pep[1:-1])-1].add(pep[0])
                    
                if all([len(wt_aa_set) == 1 for wt_aa_set in wt_aa_by_pos]):
                    self.wt_peptide = ''.join([list(aa)[0] for aa in wt_aa_by_pos])
        except:
            pass
        
        
    def t_react(self, c, K_a, n = 1, A = 1, bgrn = 0):
        return ((A-bgrn)/(1 + np.power(K_a/c, n))) + bgrn
    
    def t_react_l2_cost(self, params, r, c):
        K_a, n, A, bgrn  = params
        min_c = min(c)
#       max_c = max(c)
        
        if np.log10(K_a) < min_c:
            K_a_reg_cost = np.sqrt(self.K_a_reg)*(min_c - np.log10(K_a))
#        elif np.log10(K_a) > max_c:
#            K_a_reg_cost = np.sqrt(self.K_a_reg)*(np.log10(K_a) - max_c) 
        else:
            K_a_reg_cost = 0
            
        #Center cooperativity n at 1
        return np.concatenate([r - self.t_react(c, K_a, n = n, A = A, bgrn = bgrn), [np.sqrt(self.n_reg)*(n-1), np.sqrt(self.A_reg)*(1-A), K_a_reg_cost, np.sqrt(self.bgrn_reg)*(bgrn)  ]])
    
    def fit_indiv_t_react_curve(self, reactivity_curve):
        concs = sorted(reactivity_curve.keys())
        if len(concs) == 1: #Low reactivity at high concentration
            hill_args = [np.inf, 0, reactivity_curve[concs[0]]*2, 0]
        else:
            hill_args = opt.least_squares(self.t_react_l2_cost, 
                                       np.array([1., 1., 1., 0.]), 
                                       bounds = np.array([[0, np.inf], [0, np.inf], [0, 1], [0, 1] ]).T, 
                                       args = ([reactivity_curve[c] for c in concs], concs))['x']
        if all([r < 0.3 for r in reactivity_curve.values()]) and hill_args[0] < self.log_Ka_range[0]:
            hill_args = [np.inf, 0, reactivity_curve[concs[0]]*2, 0]
        return {'K_a': hill_args[0], 'n': hill_args[1], 'A': hill_args[2], 'bgrn': hill_args[3]}
    
    def fit_mut_pep_avidity(self):
        for pep, c_react in self.peptide_reactivity_data.items():
            self.peptide_avidity_fits[pep] = self.fit_indiv_t_react_curve(c_react)
    
    def calc_obs_C(self, pep, base_pep = '', log_units = np.exp(1)):
        if pep not in self.peptide_avidity_fits:
            return None
        else:
            c_K_a = np.clip(np.log10(self.peptide_avidity_fits[pep]['K_a']), self.log_Ka_range[0], self.log_Ka_range[1])
        if base_pep in self.peptide_avidity_fits:
            base_K_a = np.clip(np.log10(self.peptide_avidity_fits[base_pep]['K_a']), self.log_Ka_range[0], self.log_Ka_range[1])
        elif self.wt_avidity is not None:
            base_K_a = np.log10(self.wt_avidity)
        else:
            return None
        # return (c_K_a - base_K_a)/np.log10(log_units)
        return c_K_a/np.log10(log_units)
    
    def return_all_obs_C(self, log_units = np.exp(1)):
        return {pep: self.calc_obs_C(pep, log_units = log_units) for pep in self.peptide_avidity_fits}


if __name__ == '__main__':

    data_folder = './'

    K_a_reg = 0.001
    log_Ka_range = [-4, 4]

    all_tcr_info = {}

    # make the reactivity curves tsv files on the fly, extracting them from avidity_curves.xlsx
    # construct the objects
    avidity_curves = pd.read_excel(os.path.join(data_folder, 'avidity_curves.xlsx')).astype(str)
    columns_dict = {
        'tcr1': [0, 1, 2, 3],
        'tcr2': [5, 6, 7, 8],
        'tcr3': [10, 11, 12, 13],
        'tcr4': [15, 16, 17, 18],
        'tcr5': [20, 21, 22, 23],
        'tcr6': [25, 26, 27, 28],
        'tcr7': [30, 31, 32, 33]
    }

    for tcr, columns in columns_dict.items():
        tcr_avidity = avidity_curves.iloc[2:, columns]
        tcr_avidity.columns = avidity_curves.iloc[1, columns]
        tcr_avidity.to_csv(os.path.join(data_folder, tcr + '_reactivity_curves.tsv'), sep='\t', index=False)
        all_tcr_info[tcr] = TCRpMHCAvidity(mut_pep_avidity_infile = os.path.join(data_folder, tcr + '_reactivity_curves.tsv'), K_a_reg = K_a_reg, log_Ka_range= log_Ka_range)
        os.remove(os.path.join(data_folder, tcr + '_reactivity_curves.tsv'))


    ## now make the csv files

    wt_pdbs = [
        '5d2n-filtered',
        'af3_tcr2',
        'af3_tcr3',
        'af3_tcr4',
        'af3_tcr5',
        'af3_tcr6',
        'af3_tcr7'
    ]
    pep_chains = ['I', 'B', 'B', 'B', 'B', 'B', 'B']
    pep_resnum_start_list = [1, 1, 1, 1, 1, 1, 1]
    tcrs = ['tcr1', 'tcr2', 'tcr3', 'tcr4', 'tcr5', 'tcr6', 'tcr7']
    wt_peptides = ['NLVPMVATV', 'NLVPMVATV', 'NLVPMVATV', 'IMDQVPFSV', 'IMDQVPFSV', 'IMDQVPFSV', 'GRLKALCQR']

    for wt_pdb, pep_chain, pep_resnum_start, tcr, wt_pep in zip(wt_pdbs, pep_chains, pep_resnum_start_list, tcrs, wt_peptides):

        columns = ['mutant', 'is_wt', '-log10(EC50)', 'wt_pdb', 'mt_pdb', 'mutant_chain']

        df = pd.DataFrame(columns=columns)

        mutant_to_ec50 = all_tcr_info[tcr].return_all_obs_C()

        mutants = []
        is_wt = []
        ec50s = []
        sequences = []
        for mutant in mutant_to_ec50:
            wt_aa = mutant[0]
            resnum = int(mutant[1:-1])
            mt_aa = mutant[-1]

            mutant_in_pdb = f'{wt_aa}{resnum - 1 + pep_resnum_start}{mt_aa}'

            mutants.append(mutant_in_pdb)
            is_wt.append(wt_aa == mt_aa)
            ec50s.append(mutant_to_ec50[mutant])
            sequences.append(f'{wt_pep[:resnum - 1]}{mt_aa}{wt_pep[resnum:]}')

        df['mutant'] = mutants
        df['is_wt'] = is_wt
        df['-log10(EC50)'] = -np.array(ec50s)

        df['wt_seq'] = np.full(len(mutants), wt_pep)

        df['sequence'] = sequences

        # make the reliability mask
        df['is_reliable'] = df['-log10(EC50)'] >= -9.0 # this is specific to our curve fitting

        df['wt_pdb'] = np.full(len(mutants), wt_pdb)

        df['mutant_chain'] = np.full(len(mutants), pep_chain)

        df.to_csv(f'mskcc_{tcr}_ec50_sat_mut_af3.csv', index=False)

