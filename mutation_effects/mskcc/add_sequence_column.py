
import pandas as pd

from constants import TCR_TO_PDB, PDB_TO_PEP_INFO

for tcr in TCR_TO_PDB:
    pdb = TCR_TO_PDB[tcr]
    resnum_start = PDB_TO_PEP_INFO[pdb][1]
    wt_seq = PDB_TO_PEP_INFO[pdb][2]
    df = pd.read_csv(f'mskcc_tcr{tcr}_ec50_sat_mut.csv')

    df['sequence'] = [wt_seq[:int(mut[1:-1])-resnum_start] + mut[-1] + wt_seq[int(mut[1:-1])-resnum_start+1:] for mut in df['mutant']]

    df.to_csv(f'mskcc_tcr{tcr}_ec50_sat_mut.csv', index=False)
