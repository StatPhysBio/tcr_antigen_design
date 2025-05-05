
import os
import numpy as np
import pandas as pd

pdb = '2bnq'
chain = 'C'

## hermes_py_000 - fixed
df = pd.read_csv('hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t')
df['pdb_for_hermes_scoring'] = [pdb for _ in range(len(df))]
df['chain_for_hermes_scoring'] = [chain for _ in range(len(df))]
df.to_csv('hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', index=None, sep='\t')


## hermes_py_050 - fixed
df = pd.read_csv('hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', sep='\t')
df['pdb_for_hermes_scoring'] = [pdb for _ in range(len(df))]
df['chain_for_hermes_scoring'] = [chain for _ in range(len(df))]
df.to_csv('hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv', index=None, sep='\t')


## hermes_py_000 - relaxed
df = pd.read_csv('hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t')
df['pdb_for_hermes_scoring'] = [pdb for _ in range(len(df))]
df['chain_for_hermes_scoring'] = [chain for _ in range(len(df))]
df.to_csv('hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', index=None, sep='\t')


## hermes_py_050 - relaxed
df = pd.read_csv('hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t')
df['pdb_for_hermes_scoring'] = [pdb for _ in range(len(df))]
df['chain_for_hermes_scoring'] = [chain for _ in range(len(df))]
df.to_csv('hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', index=None, sep='\t')


