
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

## NOTE: excluding cysteine to compare against AW in analysis, since int their analysis they exclude cysteine

BLOSUM62 = {
    'C': {'C': 9, 'S': -1, 'T': -1, 'A': 0, 'G': -3, 'P': -3, 'D': -3, 'E': -4, 'Q': -3, 'N': -3, 'H': -3, 'R': -3, 'K': -3, 'M': -1, 'I': -1, 'L': -1, 'V': -1, 'W': -2, 'Y': -2, 'F': -2},
    'S': {'C': -1, 'S': 4, 'T': 1, 'A': 1, 'G': 0, 'P': -1, 'D': 0, 'E': 0, 'Q': 0, 'N': 1, 'H': -1, 'R': -1, 'K': 0, 'M': -1, 'I': -2, 'L': -2, 'V': -2, 'W': -3, 'Y': -2, 'F': -2},
    'T': {'C': -1, 'S': 1, 'T': 5, 'A': 0, 'G': -2, 'P': -1, 'D': -1, 'E': -1, 'Q': -1, 'N': 0, 'H': -2, 'R': -1, 'K': -1, 'M': -1, 'I': -1, 'L': -1, 'V': 0, 'W': -2, 'Y': -2, 'F': -2},
    'A': {'C': 0, 'S': 1, 'T': 0, 'A': 4, 'G': 0, 'P': -1, 'D': -2, 'E': -1, 'Q': -1, 'N': -2, 'H': -2, 'R': -1, 'K': -1, 'M': -1, 'I': -1, 'L': -1, 'V': 0, 'W': -3, 'Y': -2, 'F': -2},
    'G': {'C': -3, 'S': 0, 'T': -2, 'A': 0, 'G': 6, 'P': -2, 'D': -1, 'E': -2, 'Q': -2, 'N': 0, 'H': -2, 'R': -2, 'K': -2, 'M': -3, 'I': -4, 'L': -4, 'V': -3, 'W': -2, 'Y': -3, 'F': -3},
    'P': {'C': -3, 'S': -1, 'T': -1, 'A': -1, 'G': -2, 'P': 7, 'D': -1, 'E': -1, 'Q': -1, 'N': -1, 'H': -2, 'R': -2, 'K': -1, 'M': -2, 'I': -3, 'L': -3, 'V': -2, 'W': -4, 'Y': -3, 'F': -4},
    'D': {'C': -3, 'S': 0, 'T': -1, 'A': -2, 'G': -1, 'P': -1, 'D': 6, 'E': 2, 'Q': 0, 'N': 1, 'H': -1, 'R': -2, 'K': -1, 'M': -3, 'I': -3, 'L': -4, 'V': -3, 'W': -4, 'Y': -3, 'F': -3},
    'E': {'C': -4, 'S': 0, 'T': -1, 'A': -1, 'G': -2, 'P': -1, 'D': 2, 'E': 5, 'Q': 2, 'N': 0, 'H': 0, 'R': 0, 'K': 1, 'M': -2, 'I': -3, 'L': -3, 'V': -2, 'W': -3, 'Y': -2, 'F': -3},
    'Q': {'C': -3, 'S': 0, 'T': -1, 'A': -1, 'G': -2, 'P': -1, 'D': 0, 'E': 2, 'Q': 5, 'N': 0, 'H': 0, 'R': 1, 'K': 1, 'M': 0, 'I': -3, 'L': -2, 'V': -2, 'W': -2, 'Y': -1, 'F': -3},
    'N': {'C': -3, 'S': 1, 'T': 0, 'A': -2, 'G': 0, 'P': -2, 'D': 1, 'E': 0, 'Q': 0, 'N': 6, 'H': 1, 'R': 0, 'K': 0, 'M': -2, 'I': -3, 'L': -3, 'V': -3, 'W': -4, 'Y': -2, 'F': -3},
    'H': {'C': -3, 'S': -1, 'T': -2, 'A': -2, 'G': -2, 'P': -2, 'D': -1, 'E': 0, 'Q': 0, 'N': 1, 'H': 8, 'R': 0, 'K': -1, 'M': -2, 'I': -3, 'L': -3, 'V': -3, 'W': -2, 'Y': 2, 'F': -1},
    'R': {'C': -3, 'S': -1, 'T': -1, 'A': -1, 'G': -2, 'P': -2, 'D': -2, 'E': 0, 'Q': 1, 'N': 0, 'H': 0, 'R': 5, 'K': 2, 'M': -1, 'I': -3, 'L': -2, 'V': -3, 'W': -3, 'Y': -2, 'F': -3},
    'K': {'C': -3, 'S': 0, 'T': -1, 'A': -1, 'G': -2, 'P': -1, 'D': -1, 'E': 1, 'Q': 1, 'N': 0, 'H': -1, 'R': 2, 'K': 5, 'M': -1, 'I': -3, 'L': -2, 'V': -2, 'W': -3, 'Y': -2, 'F': -3},
    'M': {'C': -1, 'S': -1, 'T': -1, 'A': -1, 'G': -3, 'P': -2, 'D': -3, 'E': -2, 'Q': 0, 'N': -2, 'H': -2, 'R': -1, 'K': -1, 'M': 5, 'I': 1, 'L': 2, 'V': 1, 'W': -1, 'Y': -1, 'F': 0},
    'I': {'C': -1, 'S': -2, 'T': -1, 'A': -1, 'G': -4, 'P': -3, 'D': -3, 'E': -3, 'Q': -3, 'N': -3, 'H': -3, 'R': -3, 'K': -3, 'M': 1, 'I': 4, 'L': 2, 'V': 3, 'W': -3, 'Y': -1, 'F': 0},
    'L': {'C': -1, 'S': -2, 'T': -1, 'A': -1, 'G': -4, 'P': -3, 'D': -4, 'E': -3, 'Q': -2, 'N': -3, 'H': -3, 'R': -2, 'K': -2, 'M': 2, 'I': 2, 'L': 4, 'V': 1, 'W': -2, 'Y': -1, 'F': 0},
    'V': {'C': -1, 'S': -2, 'T': 0, 'A': 0, 'G': -3, 'P': -2, 'D': -3, 'E': -2, 'Q': -2, 'N': -3, 'H': -3, 'R': -3, 'K': -2, 'M': 1, 'I': 3, 'L': 1, 'V': 4, 'W': -3, 'Y': -1, 'F': -1},
    'W': {'C': -2, 'S': -3, 'T': -2, 'A': -3, 'G': -2, 'P': -4, 'D': -4, 'E': -3, 'Q': -2, 'N': -4, 'H': -2, 'R': -3, 'K': -3, 'M': -1, 'I': -3, 'L': -2, 'V': -3, 'W': 11, 'Y': 2, 'F': 1},
    'Y': {'C': -2, 'S': -2, 'T': -2, 'A': -2, 'G': -3, 'P': -3, 'D': -3, 'E': -2, 'Q': -1, 'N': -2, 'H': 2, 'R': -2, 'K': -2, 'M': -1, 'I': -1, 'L': -1, 'V': -1, 'W': 2, 'Y': 7, 'F': 3},
    'F': {'C': -2, 'S': -2, 'T': -2, 'A': -2, 'G': -3, 'P': -4, 'D': -3, 'E': -3, 'Q': -3, 'N': -3, 'H': -1, 'R': -3, 'K': -3, 'M': 0, 'I': 0, 'L': 0, 'V': -1, 'W': 1, 'Y': 3, 'F': 6}
}

AMINOACIDS = 'CMFILVWYAGTSNQDEHRKP'

def make_matrix(property):
    import json
    with open(f'aminoacid_values.json') as f:
        properties = json.load(f)[property]
    matrix = np.zeros((20, 20))
    for i, aa1 in enumerate(AMINOACIDS):
        for j, aa2 in enumerate(AMINOACIDS):
            matrix[i, j] = np.abs(properties[aa1] - properties[aa2])
    return matrix

models = ['hcnn_pyrosetta_noise=0p0', 'hcnn_pyrosetta_noise=0p5']

for model in models:

    matrix_num = np.zeros((20, 20))
    matrix_denom = np.zeros((20, 20))

    all_files = os.listdir(f'{model}/')

    for file in all_files:
        if file.endswith('.csv'):
            wt_pep = file.strip('.csv').split('__')[1]
            df = pd.read_csv(f'{model}/{file}')
            for i, row in df.iterrows():
                aa_wt = wt_pep[i]
                for j, aa_mt in enumerate(AMINOACIDS):
                    matrix_num[AMINOACIDS.index(aa_wt), AMINOACIDS.index(aa_mt)] += np.log(row[aa_mt]) - np.log(row[aa_wt])
                    matrix_denom[AMINOACIDS.index(aa_wt), AMINOACIDS.index(aa_mt)] += 1

    matrix = matrix_num / matrix_denom

    matrix_df = pd.DataFrame(matrix, index=list(AMINOACIDS), columns=list(AMINOACIDS))
    matrix_df.to_csv(f'{model}_marginal_matrix.csv')

    # excluding cysteines!!
    matrix_df.loc['C', :] = np.nan
    matrix_df.loc[:, 'C'] = np.nan

    plt.imshow(matrix, cmap='viridis', vmin=np.min(matrix), vmax=np.max(matrix))
    plt.xticks(np.arange(20), list(AMINOACIDS))
    plt.yticks(np.arange(20), list(AMINOACIDS))
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(f'{model}_marginal_matrix.png')
    plt.close()

    properties = ['V17', 'V3', 'V14', 'V44', 'V50', 'V49', 'V19', 'V22', 'V25']

    correlations = []
    for property in properties:
        matrix = make_matrix(property)
        for i in range(20):
            matrix[i, i] = np.nan
        matrix_hermes = matrix_df.values
        for i in range(20):
            matrix_hermes[i, i] = np.nan
        
        matrix_flat = matrix.flatten()

        matrix_hermes_flat = matrix_hermes.flatten()

        mask = ~np.isnan(matrix_flat) & ~np.isnan(matrix_hermes_flat)
        matrix_flat = matrix_flat[mask]
        matrix_hermes_flat = matrix_hermes_flat[mask]
        
        from scipy.stats import spearmanr
        sr, sr_pval = spearmanr(-matrix_flat, matrix_hermes_flat)
        print(f'{model}, {property}, {sr:.2f}')
        correlations.append(sr)
    print()


    xticks = [0, 1, 2, 4, 5, 6, 8, 9, 10]
    plt.bar(xticks, correlations)
    plt.xticks([0, 1, 2, 4, 5, 6, 8, 9, 10], properties)
    plt.ylabel('Spearman correlation')
    plt.tight_layout()
    plt.savefig(f'{model}_marginal_matrix_correlations.png')
    plt.close()




    blosum_matrix = np.zeros((20, 20))
    for i, aa1 in enumerate(AMINOACIDS):
        for j, aa2 in enumerate(AMINOACIDS):
            blosum_matrix[i, j] = BLOSUM62[aa1][aa2]
    # exclude diagonal, set to nan
    for i in range(20):
        blosum_matrix[i, i] = np.nan
    
    blosum_matrix_flat = blosum_matrix.flatten()
    blosum_matrix_flat = blosum_matrix_flat[~np.isnan(blosum_matrix_flat)]


    correlations = []
    for property in properties:
        matrix = make_matrix(property)
        for i in range(20):
            matrix[i, i] = np.nan
        matrix_hermes = matrix_df.values
        for i in range(20):
            matrix_hermes[i, i] = np.nan
        
        matrix_flat = matrix.flatten()
        matrix_flat = matrix_flat[~np.isnan(matrix_flat)]

        from scipy.stats import spearmanr
        sr, sr_pval = spearmanr(-matrix_flat, blosum_matrix_flat)
        print(f'{model}, {property}, {sr:.2f}')
        correlations.append(sr)
    print()

    xticks = [0, 1, 2, 4, 5, 6, 8, 9, 10]
    plt.bar(xticks, correlations)
    plt.xticks([0, 1, 2, 4, 5, 6, 8, 9, 10], properties)
    plt.ylabel('Spearman correlation')
    plt.tight_layout()
    plt.savefig(f'{model}_blosum_correlations.png')
    plt.close()


    aw_matrix = pd.read_csv('../AW_d.csv', sep='\t', index_col=0)

    values_aw = []
    values_hermes = []

    for aa1 in AMINOACIDS:
        for aa2 in AMINOACIDS:
            if aa1 != aa2 and aa1 in aw_matrix.index and aa2 in aw_matrix.columns:

                value_aw = aw_matrix.loc[aa1, aa2]
                value_hermes = matrix_df.loc[aa1, aa2]

                if not np.isnan(value_aw) and not np.isnan(value_hermes):
                    values_aw.append(aw_matrix.loc[aa1, aa2])
                    values_hermes.append(matrix_df.loc[aa1, aa2])

    print(spearmanr(values_aw, values_hermes))
    print()
    print()



