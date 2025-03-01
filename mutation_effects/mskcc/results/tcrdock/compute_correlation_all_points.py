
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

for tcr in range(1, 8):
    df = pd.read_csv(f'mskcc_tcr{tcr}_ec50_sat_mut_af3_w_pae.tsv', sep='\t')
    print(tcr, spearmanr(df['-log10(EC50)'], -df['pmhc_tcr_pae']))
