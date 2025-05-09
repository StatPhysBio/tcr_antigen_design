
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

scores_true = pd.read_csv('hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae.tsv', sep='\t')['pnE'].values
scores_remake = pd.read_csv('hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv', sep='\t')['pnE'].values

plt.scatter(scores_true, scores_remake)
plt.xlabel('True')
plt.ylabel('Remake')

xlim = plt.gca().get_xlim()
ylim = plt.gca().get_ylim()

lim = (min(xlim[0], ylim[0]), max(xlim[1], ylim[1]))

plt.xlim(lim)
plt.ylim(lim)

plt.plot(lim, lim, color='black')

plt.title(f'Pr: %.2f' % (pearsonr(scores_true, scores_remake)[0]))

plt.savefig('test_with_max.png')
plt.close()


