import os

from constants import WT_SEQ

os.system(f'python ../../src/tcrdock_and_blosum_plots.py --wt_seq {WT_SEQ}')
