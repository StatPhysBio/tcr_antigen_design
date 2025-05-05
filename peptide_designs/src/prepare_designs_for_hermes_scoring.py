
'''

Score every peptide on every structure and with every protocol

'''

import os
import numpy as np
import pandas as pd

    #             --model_version $model_version \
    #             --pdbdir ./pdbs \
    #             --output_csv_filepath $output_dir$model_version'/nyeso_peptide_kd_closest-'$model_version'-use_mt_structure=0.csv' \
    #             --csv_filepath nyeso_peptide_kd_closest.csv \
    #             --pdb_column pdb \
    #             --chain_column chain \
    #             --peptide_column sequence \
    #             --peptide_resnum_start 1


## nyeso
df_in = pd.read_excel('../All-designs.xlsx', sheet_name='NYES0-designs')

sequences = df_in['sequence'].values
sequences = list(filter(lambda x: isinstance(x, str), sequences))

pdbs = ['2bnr', '2bnq']
chains = ['C', 'C']

dict_out = {
    'sequence': [],
    'pdb': [],
    'chain': []
}

for seq in sequences:
    for pdb, chain in zip(pdbs, chains):
        dict_out['sequence'].append(seq)
        dict_out['pdb'].append(pdb)
        dict_out['chain'].append(chain)

df_out = pd.DataFrame(dict_out)

df_out.to_csv('../nyeso/selected_designs_for_hermes_scoring.csv', index=None)
        

## ebv
df_in = pd.read_excel('../All-designs.xlsx', sheet_name='EBV-designs')

sequences = df_in['sequence'].values
sequences = list(filter(lambda x: isinstance(x, str), sequences))

pdbs = ['3mv7.pdb.human.MH1.B-35.A.C.DE', '4prp.pdb.human.MH1.B-35.A.C.DE']
chains = ['B', 'B']

dict_out = {
    'sequence': [],
    'pdb': [],
    'chain': []
}

for seq in sequences:
    for pdb, chain in zip(pdbs, chains):
        dict_out['sequence'].append(seq)
        dict_out['pdb'].append(pdb)
        dict_out['chain'].append(chain)

df_out = pd.DataFrame(dict_out)

df_out.to_csv('../ebv/selected_designs_for_hermes_scoring.csv', index=None)



## mage
df_in = pd.read_excel('../All-designs.xlsx', sheet_name='MAGE-designs')

sequences = df_in['sequence'].values
sequences = list(filter(lambda x: isinstance(x, str), sequences))

pdbs = ['5brz.pdb.human.MH1.A-01.A.C.DE', '5bs0.pdb.human.MH1.A-01.A.C.DE']
chains = ['B', 'B']

dict_out = {
    'sequence': [],
    'pdb': [],
    'chain': []
}

for seq in sequences:
    for pdb, chain in zip(pdbs, chains):
        dict_out['sequence'].append(seq)
        dict_out['pdb'].append(pdb)
        dict_out['chain'].append(chain)

df_out = pd.DataFrame(dict_out)

df_out.to_csv('../magea3_and_titin/selected_designs_for_hermes_scoring.csv', index=None)


