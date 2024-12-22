#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 22:34:26 2020

@author: zacharysethna

Edited Sat Dec 21 2024

@author: Armita Nourmohammad

"""




import numpy as np
import scipy.optimize as opt
from itertools import chain, combinations
import time
import json
import os
from scipy.stats import linregress

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
        
        
#%%
class InferEpitopeDist(object):
    def __init__(self, blosum62_json_file = None, blosum62_reg = 0):
        self.amino_acids = 'ACDEFGHIKLMNPQRSTVWY'
        self.amino_acid_dict = {aa: i for i, aa in enumerate(self.amino_acids)}

        self.blosum62_reg = blosum62_reg
        self.embed_bound = 10
        self.n_dim = 2
        self.log_units = np.exp(1)
        #self.torus_topology = [False for _ in range(self.n_dim)]
        self.tcrs = {}
        self.reactivity_cutoff = 0.1
        self.inferred_epitope_dist_models = []
        
        if blosum62_json_file is None:
            blosum62_json_file = os.path.join(os.path.dirname(__file__), 'blosum62.json')
        with open(blosum62_json_file, 'r') as infile:
            self.blosum62_dict = json.load(infile)
        blosum62_mat = np.zeros((len(self.amino_acids), len(self.amino_acids)))
        for sub, val in self.blosum62_dict.items():
            try:
                blosum62_mat[self.amino_acid_dict[sub[0]], self.amino_acid_dict[sub[-1]]] = val
                blosum62_mat[self.amino_acid_dict[sub[-1]], self.amino_acid_dict[sub[0]]] = val
            except KeyError:
                pass
        self.blosum62_mat = blosum62_mat
        self.blosum62_reg_mask = np.logical_not(np.eye(len(self.amino_acids)))
        self.blosum62_var = np.sum((self.blosum62_mat[self.blosum62_reg_mask] - np.mean(self.blosum62_mat[self.blosum62_reg_mask]))**2)            

    def ham1pep2p_inds(self, ham1pep):
        '''Convert hamming dist 1 peptide representation p to index representation'''
        return [int(x) for x in [int(ham1pep[1:-1])-1, self.amino_acid_dict[ham1pep[0]], self.amino_acid_dict[ham1pep[-1]]]]
    
    def ham1tcr_dist(self, p_inds, d_i, M_ab):
        return d_i[p_inds[0]]*M_ab[p_inds[1], p_inds[2]]
    
    def calc_energy_cost(self, d_i, M_ab, p_dist_arr):
        '''p_dist_arr is n x 4 array of i, aa_1, aa_2, p_dist for each peptide'''
        
    #    E = np.sum([(p_dists[p_inds] - ham1tcr_dist(p_inds, d_i, M_ab))**2 for p_inds in p_dists])/len(p_dists)
    #    return E
        E = np.sum([(c_p_dist[3] - d_i[c_p_dist[0]]*M_ab[c_p_dist[1], c_p_dist[2]])**2 for c_p_dist in p_dist_arr])/len(p_dist_arr)
        if self.blosum62_reg == 0:
            return E
        else:
            slope_fit, intercept_fit, r_value_fit, p_value_fit, std_err_fit = linregress(np.array([M_ab[self.blosum62_reg_mask], self.blosum62_mat[self.blosum62_reg_mask]]).T)
            return E + self.blosum62_reg*(1-r_value_fit**2)*self.blosum62_var

    def calc_centered_d_i(self, variance):
        d_i = np.exp(-(np.arange(9)-4)**2/(2*variance))
        d_i = d_i/np.sum(d_i)
        return d_i
        
    def calc_M_ab_euclid_coords(self, euclid_coords):
        n_aa = euclid_coords.shape[0]
        M_ab = np.zeros((n_aa, n_aa))
        for i in range(euclid_coords.shape[0]):
            for j in range(i+1, euclid_coords.shape[0]):
                M_ab[i, j] = np.linalg.norm(euclid_coords[i] - euclid_coords[j])
                M_ab[j, i] = M_ab[i, j]      
        return M_ab

    def process_euclid_arg_input(self, x):
        d_i = x[:9]
        euclid_coords = x[9:].reshape(len(self.amino_acids), self.n_dim)
#        d_i = self.calc_centered_d_i(x[0])
#        euclid_coords = x[1:].reshape(len(self.amino_acids), self.n_dim)
        return d_i, euclid_coords
    
    def euclid_tcr_xrx_cost_func(self, x, p_dist_arr):
        d_i, euclid_coords = self.process_euclid_arg_input(x)
        M_ab = self.calc_M_ab_euclid_coords(euclid_coords)
        return self.calc_energy_cost(d_i, M_ab, p_dist_arr)
    
    def make_p_dist_arr_from_tcr(self, tcr):
        #tcr is of TCRpMHCAvidity class
        p_dist_arr = []
        for pep, c_fit in tcr.peptide_avidity_fits.items():
            #Exclude wt peptides for speed (constant cost offset should be irrelevant either way)
            if pep[0] == pep[-1]: continue 
            c_p_inds = self.ham1pep2p_inds(pep)
            c_dist = tcr.calc_obs_C(pep, log_units = self.log_units)
            p_dist_arr.append(sum([c_p_inds, [c_dist]], []))
        return p_dist_arr
    
    def make_p_dist_arr_all_combos_from_tcr(self, tcr):
        #tcr is of TCRpMHCAvidity class
        p_dist_arr = []
        c_peptides = list(tcr.peptide_avidity_fits.keys())
        pep_mut_pos = sorted(set([int(pep[1:-1]) for pep in c_peptides]))
        peptides_by_positions = {c_pos: [] for c_pos in pep_mut_pos}
        for pep in tcr.peptide_avidity_fits.keys():
            c_pos = int(pep[1:-1])
            peptides_by_positions[c_pos].append(pep)
        for c_pos, c_peps in peptides_by_positions.items():
            for i, pepA in enumerate(c_peps):
                for pepB in c_peps[i+1:]:
                    if (tcr.peptide_avidity_fits[pepA]['K_a'] > self.reactivity_cutoff) and (tcr.peptide_avidity_fits[pepB]['K_a'] > self.reactivity_cutoff): continue
                    c_p_inds = self.ham1pep2p_inds(pepA[-1] + pepB[1:])
                    c_dist = tcr.calc_obs_C(pepB, base_pep = pepA, log_units = self.log_units)
                    p_dist_arr.append(sum([c_p_inds, [c_dist]], []))
        return p_dist_arr

    def infer_euclid_coord_model_from_p_dist_arr(self, p_dist_arr, **kwarg):
        if 'seed' in kwarg: np.random.seed(kwarg['seed'])
        c_bounds = np.zeros((9+len(self.amino_acids)*self.n_dim, 2))
        c_bounds[:9, 1] = 100/self.embed_bound
        c_bounds[9:, 1] = self.embed_bound
        c_bounds[9:, 0] = -self.embed_bound
        
#        c_bounds = np.zeros((1+len(self.amino_acids)*self.n_dim, 2))
#        c_bounds[0, 1] = 20
#        c_bounds[1:, 1] = self.embed_bound
#        c_bounds[1:, 0] = -self.embed_bound
        
        start_time = time.time()
        model_fit = opt.dual_annealing(self.euclid_tcr_xrx_cost_func, bounds = c_bounds, args = (p_dist_arr,), maxfun = int(1e6))
        d_i, euclid_coords = self.process_euclid_arg_input(model_fit.x)
        print('Completed model inference in %.2f seconds'%(time.time() - start_time))
        return {'d_i': d_i, 'euclid_coords': euclid_coords, 'M_ab': self.calc_M_ab_euclid_coords(euclid_coords)}
    
    def infer_euclid_coord_model_from_tcrs(self, tcrs, **kwarg):
        if 'blosum62_reg' in kwarg: self.blosum62_reg = kwarg['blosum62_reg']
        if type(tcrs) is not list: tcrs = [tcrs]
        if 'use_all_pep_combos' in kwarg:
            if kwarg['use_all_pep_combos']:
                c_p_dist_arr = sum([self.make_p_dist_arr_all_combos_from_tcr(tcr) for tcr in tcrs], [])
            else:
                c_p_dist_arr = sum([self.make_p_dist_arr_from_tcr(tcr) for tcr in tcrs], [])
        else:
            c_p_dist_arr = sum([self.make_p_dist_arr_from_tcr(tcr) for tcr in tcrs], [])
        print('Starting model inference...')
        c_model = self.infer_euclid_coord_model_from_p_dist_arr(c_p_dist_arr, **kwarg)
        c_model['blosum62_reg'] = self.blosum62_reg
        c_model['log_units'] = self.log_units
        self.inferred_epitope_dist_models.append(c_model)
        return c_model
    
    def calc_M_ab_bl62(self, poly_coords):
        M_ab = sum([x*(self.blosum62_mat**i) for i, x in enumerate(poly_coords)])*(np.ones(self.blosum62_mat.shape) - np.eye(len(self.amino_acids)))
        return M_ab
    
    def process_bl62_arg_input(self, x):
        d_i = x[:9]
        poly_coords = x[9:]
        
#        d_i = self.calc_centered_d_i(x[0])
#        poly_coords = x[1:]
        return d_i, poly_coords
    
    def bl62_tcr_xrx_cost_func(self, x, p_dist_arr):
        d_i, poly_coords = self.process_bl62_arg_input(x)
        M_ab = self.calc_M_ab_bl62(poly_coords)
        return self.calc_energy_cost(d_i, M_ab, p_dist_arr)
    
    def infer_bl62_model_from_p_dist_arr(self, p_dist_arr, **kwarg):
        if 'seed' in kwarg: np.random.seed(kwarg['seed'])
        if 'poly_fit_dim' in kwarg: 
            poly_fit_dim = kwarg['poly_fit_dim']
        else:
            poly_fit_dim = 2
        c_bounds = np.zeros((9+poly_fit_dim, 2))
        c_bounds[:9, 0] = 1e-3
        c_bounds[:9, 1] = 10
        c_bounds[9:, 0] = -200
        c_bounds[9:, 1] = 200
        
#        c_bounds = np.zeros((1+poly_fit_dim, 2))
#        c_bounds[0, 1] = 20
#        c_bounds[1:, 0] = -20
#        c_bounds[1:, 1] = 20
        
        start_time = time.time()
        model_fit = opt.dual_annealing(self.bl62_tcr_xrx_cost_func, bounds = c_bounds, args = (p_dist_arr,), maxfun = int(1e6))
        d_i, poly_coords = self.process_bl62_arg_input(model_fit.x)
        print('Completed model inference in %.2f seconds'%(time.time() - start_time))
        return {'d_i': d_i, 'poly_coords': poly_coords, 'M_ab': self.calc_M_ab_bl62(poly_coords), "model_fit": model_fit.x}
    
    def infer_bl62_model_from_tcrs(self, tcrs, **kwarg):
        if type(tcrs) is not list: tcrs = [tcrs]
        if 'use_all_pep_combos' in kwarg:
            if kwarg['use_all_pep_combos']:
                c_p_dist_arr = sum([self.make_p_dist_arr_all_combos_from_tcr(tcr) for tcr in tcrs], [])
            else:
                c_p_dist_arr = sum([self.make_p_dist_arr_from_tcr(tcr) for tcr in tcrs], [])
        else:
            c_p_dist_arr = sum([self.make_p_dist_arr_from_tcr(tcr) for tcr in tcrs], [])
        print('Starting model inference...')
        c_model = self.infer_bl62_model_from_p_dist_arr(c_p_dist_arr, **kwarg)
        c_model['log_units'] = self.log_units
        self.inferred_epitope_dist_models.append(c_model)
        return c_model