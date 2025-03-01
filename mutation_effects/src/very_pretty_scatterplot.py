
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, spearmanr

import argparse

import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams['figure.facecolor'] = 'white'

IND_TO_AA = 'CSTAGPDEQNHRKMILVWYF'
AA_TO_IND = {aa: i for i, aa in enumerate(IND_TO_AA)}

def system_name_to_wt_seq(system_name_in_csv_file):
    if 'nyeso' in system_name_in_csv_file:
        return 'SLLMWITQC'
    elif 'tax' in system_name_in_csv_file:
        return 'LLFGYPVYV'
    elif 'hsiue' in system_name_in_csv_file:
        return 'HMTEVVRHC'
    elif 'mskcc' in system_name_in_csv_file:
        if 'tcr1' in system_name_in_csv_file or 'tcr2' in system_name_in_csv_file or 'tcr3' in system_name_in_csv_file:
            return 'NLVPMVATV'
        elif 'tcr4' in system_name_in_csv_file or 'tcr5' in system_name_in_csv_file or 'tcr6' in system_name_in_csv_file:
            return 'IMDQVPFSV'
        elif 'tcr7' in system_name_in_csv_file:
            return 'GRLKALCQR'
    else:
        raise ValueError(f'Unknown system: {system_name_in_csv_file}')

PDB_TO_COLOR = {
    '2bnq': '#4682B4',
    '2bnr': '#FF6F61',
    '1ao7': '#4682B4',
    '1qsf': '#FF6F61',
}


SYSTEM_TO_WT_TO_SHOW_COLUMNS = {
    'nyeso': None,
    'tax': None,
    'mart': None,
    'hsiue_et_al': 'is_wt',
    'mskcc': 'is_wt'
}

SYSTEM_TO_TARGET_COLUMN = {
    'nyeso': '-log10(Kd)',
    'tax': '-log10(Kd)',
    'mart': '-log10(Kd)',
    'hsiue_et_al': 'IFN_gamma (pg/ml)',
    'mskcc': '-log10(EC50)'
}

# SYSTEM_TO_PRETTY_TARGET = {
#     'nyeso': '$-\\text{log}_{10}(K_d)$',
#     'tax': '$-\\text{log}_{10}(K_d)$',
#     'mart': '$-\\text{log}_{10}(K_d)$',
#     'hsiue_et_al': 'IFN_gamma (pg/ml)',
#     'mskcc': '$- \Delta \\text{log(EC50)}$'
# }

SYSTEM_TO_PRETTY_TARGET = {
    'nyeso': '$-\\text{log}_{10}(K_d)$',
    'tax': '$-\\text{log}_{10}(K_d)$',
    'mart': '$-\\text{log}_{10}(K_d)$',
    'hsiue_et_al': 'IFN_gamma (pg/ml)',
    'mskcc': '$-\\text{log}_{10}(\\text{EC}_{50})$'
}

SYSTEM_TO_MARKER_SIZE = {
    'nyeso': 60,
    'tax': 60,
    'mart': 60,
    'hsiue_et_al': 30,
    'mskcc': 30
}

SYSTEM_TO_ALPHA = {
    'nyeso': 0.6,
    'tax': 0.6,
    'mart': 0.6,
    'hsiue_et_al': 0.6,
    'mskcc': 0.6
}

MOTIF_TO_MARKER = {
    'GIG': 'd',
    'DRG': 'x'
}

SYSTEM_CSV_TO_TITLE = {
    'nyeso_peptide_kd_closest': 'NY-ESO-1',
    'tax_peptide_kd_closest': 'Tax',
    'mart_peptide_kd_closest': 'MART-1 (right motifs)',
    'mart_peptide_kd_3qdg': 'MART-1 (3qdg, GIG motif)',
    'mart_peptide_kd_6am5': 'MART-1 (6am5, GIG motif)',
    'mart_peptide_kd_6amu': 'MART-1 (6amu, DRG motif)',
    'hsiue_et_al_H2_sat_mut': 'H2-scDb on p53$^\\text{R175H}$',
    'mskcc_tcr1_ec50_sat_mut_af3': 'TCR1 on NLVPMVATV',
    'mskcc_tcr2_ec50_sat_mut_af3': 'TCR2 on NLVPMVATV (AF3)',
    'mskcc_tcr3_ec50_sat_mut_af3': 'TCR3 on NLVPMVATV (AF3)',
    'mskcc_tcr4_ec50_sat_mut_af3': 'TCR4 on IMDQVPFSV (AF3)',
    'mskcc_tcr5_ec50_sat_mut_af3': 'TCR5 on IMDQVPFSV (AF3)',
    'mskcc_tcr6_ec50_sat_mut_af3': 'TCR6 on IMDQVPFSV (AF3)',
    'mskcc_tcr7_ec50_sat_mut_af3': 'TCR7 on GRLKALCQR (AF3)',
}

'''
Scatterplot for chosen HERMES model
'''



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_version', type=str, required=True)
    parser.add_argument('--with_relaxation', type=int, choices=[0, 1], required=True)
    parser.add_argument('--system', type=str, required=True)
    parser.add_argument('--system_name_in_csv_file', type=str, required=True)

    args = parser.parse_args()

    target_column = SYSTEM_TO_TARGET_COLUMN[args.system]
    prediction_column = 'pnE'
    title = SYSTEM_CSV_TO_TITLE[args.system_name_in_csv_file]
    xlabel = SYSTEM_TO_PRETTY_TARGET[args.system]
    # ylabel = r'$E(\sigma_p ; \text{TCR}, \text{MHC})$'
    if args.with_relaxation:
        ylabel = 'HERMES-relaxed prediction'
    else:
        ylabel = 'HERMES-fixed prediction'
    color = 'tab:purple'
    s = SYSTEM_TO_MARKER_SIZE[args.system]
    alpha = SYSTEM_TO_ALPHA[args.system]
    
    if args.model_version == 'hermes_py_000' and args.with_relaxation:
        model_pretty_name = 'HERMES-relaxed 0.00'
    elif args.model_version == 'hermes_py_000' and not args.with_relaxation:
        model_pretty_name = 'HERMES-fixed 0.00'
    elif args.model_version == 'hermes_py_050' and args.with_relaxation:
        model_pretty_name = 'HERMES-relaxed 0.50'
    elif args.model_version == 'hermes_py_050' and not args.with_relaxation:
        model_pretty_name = 'HERMES-fixed 0.50'
    else:
        raise ValueError(f'Unknown model version: {args.model_version}')

    with_relaxation_str = '_with_relaxation' if args.with_relaxation else ''

    df = pd.read_csv(f'/{args.system}/results/{args.model_version}/{args.system_name_in_csv_file}{with_relaxation_str}-{args.model_version}-use_mt_structure=0.csv')

    # exlude nans, keep only one wildtype measurement
    is_not_nan_mask = np.logical_and(~np.isnan(df[target_column]), ~np.isnan(df[prediction_column]))
    
    if 'is_wt' in df.columns:
        is_wt_mask = df['is_wt'].astype(bool)
        first_wt_mask = np.zeros(df.shape[0], dtype=bool)
        first_wt_mask[np.arange(df.shape[0])[is_wt_mask][0]] = True
        mask = np.logical_or(np.logical_and(is_not_nan_mask, ~is_wt_mask), first_wt_mask)
    else:
        print('Warning: no "is_wt" column in csv file, assuming there are no wildtype duplicates.', file=sys.stderr)
        mask = is_not_nan_mask

    df = df[mask]

    predictions = df[prediction_column]
    targets = df[target_column]

    if 'motif' in df.columns:
        markers = np.array([MOTIF_TO_MARKER[motif] for motif in df['motif'].values])
    else:
        markers = np.array(['o'] * df.shape[0])
    unique_markers = np.unique(markers)

    plt.figure(figsize=(3.65, 3.65))
    
    if SYSTEM_TO_WT_TO_SHOW_COLUMNS[args.system] is not None:
        wt_values = df[df[SYSTEM_TO_WT_TO_SHOW_COLUMNS[args.system]].astype(bool)][[target_column, prediction_column]].values
        for wt_value_target, wt_value_pred in wt_values:
            plt.axvline(x=wt_value_target, c='k', alpha=1.0, ls='--', lw=0.85)
            plt.axhline(y=wt_value_pred, c='k', alpha=1.0, ls='--', lw=0.85)

    if args.system == 'mskcc':

        wt_target_value = wt_values[0][0] # there's only one anyways
        print(f'wt_target_value = {wt_target_value}')

        is_reliable_mask = df['is_reliable'].values

        targets_rel = targets[is_reliable_mask]
        predictions_rel = predictions[is_reliable_mask]

        target_non_rel = targets[~is_reliable_mask]
        prediction_non_rel = predictions[~is_reliable_mask]

        pr, pr_pval = pearsonr(targets_rel, predictions_rel)
        sr, sr_pval = spearmanr(targets_rel, predictions_rel)
        num = targets_rel.shape[0]

        from sklearn.metrics import roc_auc_score
        # auroc_all = roc_auc_score(targets > wt_target_value, predictions)
        # auroc_rel = roc_auc_score(targets_rel > wt_target_value, predictions_rel)
        try:
            auroc = roc_auc_score(df['is_reliable'].values, predictions)
        except ValueError:
            auroc = np.nan

        plt.scatter(targets_rel, predictions_rel, c=color, edgecolors="none", marker=markers[0], alpha=alpha, s=s)
        plt.scatter(target_non_rel, prediction_non_rel, c='tab:grey', edgecolors="none", marker=markers[0], alpha=alpha, s=s)

        plt.text(0.03, 0.03, f'Spearman r = {sr:.2f}\np-val = {sr_pval:.1e}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='left')
        # plt.text(0.03, 0.03, f'Sr w\out grey = {sr:.2f} ({sr_pval:.1e})\nAUROC w\out grey = {auroc_rel:.2f}\nAUROC w\ grey = {auroc_all:.2f}', transform=plt.gca().transAxes, fontsize=11, verticalalignment='bottom', horizontalalignment='left')
        # plt.text(0.03, 0.03, f'Sr = {sr:.2f} (pv {sr_pval:.1e})\nAUROC = {auroc:.2f}', transform=plt.gca().transAxes, fontsize=11, verticalalignment='bottom', horizontalalignment='left')

    elif args.system in {'nyeso', 'tax'}:

        pr, pr_pval = pearsonr(targets, predictions)
        sr, sr_pval = spearmanr(targets, predictions)
        num = targets.shape[0]

        pdbs = df['pdb'].values
        colors = [PDB_TO_COLOR[pdb] for pdb in pdbs]

        plt.scatter(targets, predictions, c=colors, edgecolors="none", alpha=alpha, s=s)

        # put text of correlation coefficient
        plt.text(0.03, 0.03, f'Spearman r = {sr:.2f}\np-val = {sr_pval:.1e}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='left')
    
    else:

        pr, pr_pval = pearsonr(targets, predictions)
        sr, sr_pval = spearmanr(targets, predictions)
        num = targets.shape[0]

        unique_markers = np.unique(markers)
        for marker in unique_markers:
            mask = markers == marker
            plt.scatter(targets[mask], predictions[mask], c=color, edgecolors="none", marker=marker, alpha=alpha, s=s)

        # put text of correlation coefficient
        plt.text(0.03, 0.03, f'Spearman r = {sr:.2f}\np-val = {sr_pval:.1e}', transform=plt.gca().transAxes, fontsize=12, verticalalignment='bottom', horizontalalignment='left')

    
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)

    plt.xlabel(xlabel, fontsize=14)
    plt.ylabel(ylabel, fontsize=14)

    plt.title(title, fontsize=15)
    # plt.title(title+'\n'+f'Spearman r = {sr:.2f} (pv {sr_pval:.1e})', fontsize=12)

    # plt.legend(loc='upper left')

    plt.tight_layout()
    outfile = f'/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-scatterplot.png'
    plt.savefig(outfile, bbox_inches="tight", dpi=300, transparent=False)
    plt.savefig(outfile.replace('.png', '.svg'), format="svg", bbox_inches="tight", dpi=300, transparent=False)
    plt.savefig(outfile.replace('.png', '.pdf'), format="pdf", bbox_inches="tight", dpi=300, transparent=False)
    plt.close()

    if len(unique_markers) > 1:
        # make legend in standlone figure
        plt.figure(figsize=(1.5, 1.5))
        for motif in MOTIF_TO_MARKER:
            marker = MOTIF_TO_MARKER[motif]
            plt.scatter([], [], c=color, marker=marker, alpha=alpha, s=s, label=motif)
        plt.legend(loc='center')
        plt.axis('off')
        plt.tight_layout()
        outfile = f'/{args.system}/plots/motif_legend.png'
        plt.savefig(outfile, bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(outfile.replace('.png', '.svg'), format="svg", bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(outfile.replace('.png', '.pdf'), format="pdf", bbox_inches="tight", dpi=300, transparent=False)
        plt.close()
    

    # now also make two heatmaps! one with the predictions, one with the targets

    MUTANT_SPLIT_SYMBOL = '|'

    def compute_mutant(wt_seq, mt_seq):
        assert len(wt_seq) == len(mt_seq), f'Lengths of sequences do not match: {wt_seq} vs {mt_seq}'
        mutants = []
        for i in range(len(wt_seq)):
            if wt_seq[i] != mt_seq[i]:
                mutants.append(f'{wt_seq[i]}{i+1}{mt_seq[i]}')
        return MUTANT_SPLIT_SYMBOL.join(mutants), len(mutants)

    def make_colorbar(data, label, filename, cmap, figsize, orientation='horizontal', limits=None):

        import matplotlib.pyplot as plt
        import matplotlib.colors as mcolors
        import matplotlib.colorbar as cbar
        import numpy as np

        # Define a colormap and normalization based on the data
        if limits is not None:
            norm = mcolors.Normalize(vmin=limits[0], vmax=limits[1])
        else:
            norm = mcolors.Normalize(vmin=np.min(data), vmax=np.max(data))

        # Standalone colorbar figure
        plt.figure(figsize=figsize)
        ax2 = plt.gca()
        # ax2.set_visible(False)  # Hide the axis
        colorbar = cbar.ColorbarBase(ax2, cmap=cmap, norm=norm, orientation=orientation)
        # colorbar.set_label(label)

        # linspaced = np.linspace(np.min(data), np.max(data), 6)[1:-1]
        # int_rounded = np.round(linspaced).astype(int)
        # colorbar.set_ticks(int_rounded)

        # # get xtick values
        # ticks = colorbar.get_ticks()
        # new_ticks = []
        # for tick in ticks:
        #     if int(tick * 10) % 10 == 0:
        #         print(tick)
        #         new_ticks.append(int(tick))
        # # set new ticks
        # colorbar.set_ticks(new_ticks)

        # set fontsize
        for t in colorbar.ax.yaxis.get_ticklabels():
            t.set_fontsize(14)

        # set fontsize
        for t in colorbar.ax.xaxis.get_ticklabels():
            t.set_fontsize(14)

        plt.tight_layout()
        plt.savefig(filename, bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(filename.replace('.png', '.svg'), format="svg", bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(filename.replace('.png', '.pdf'), format="pdf", bbox_inches="tight", dpi=300, transparent=False)
        plt.close()


    if args.system in {'nyeso', 'tax', 'mart'}:
        # make the heatmaps vertically, with the mutations on the y-axis.

        # wt_seq = system_name_to_wt_seq(args.system_name_in_csv_file)

        # sequences = df['sequence'].values

        # mutants, num_mutants = zip(*[compute_mutant(wt_seq, seq) for seq in sequences])

        # mutants = np.array(mutants)
        # num_mutants = np.array(num_mutants)

        # mutants[num_mutants == 0] = 'WT'

        # # sorting indices based on number of mutations
        # sorted_indices = np.argsort(num_mutants)

        # mutants = mutants[sorted_indices]

        # target_values = df[target_column].values[sorted_indices]
        # predicted_values = df[prediction_column].values[sorted_indices]

        mutants = df['mutants'].values
        resnums = [int(mut[1]) for mut in mutants]
        pdbs = df['pdb'].values
        
        # sort first by resnum of first mutation, then by pdb
        sorted_indices = np.lexsort((resnums, pdbs))

        target_values = df[target_column].values[sorted_indices]
        predicted_values = df[prediction_column].values[sorted_indices]
        mutants = mutants[sorted_indices]
        pdbs = pdbs[sorted_indices]
        colors = [PDB_TO_COLOR[pdb] for pdb in pdbs]

        new_mutants = []
        for mut, pdb in zip(mutants, pdbs):
            if len(mut.split(MUTANT_SPLIT_SYMBOL)) == 1 and mut[0] == mut[-1]:
                new_mutants.append(f'WT: {pdb.upper()}')
            else:
                new_mutants.append(mut)
        mutants = np.array(new_mutants)

        # make the heatmap
        figsize = (len(mutants)*0.45, 2.0)
        fontsize = 15
        rotation = 65

        plt.figure(figsize=figsize)
        heatmap = target_values[:, np.newaxis].T
        plt.imshow(heatmap, aspect='auto', cmap='viridis', alpha=1)
        plt.yticks([])
        plt.xticks(range(len(mutants)), mutants, fontsize=fontsize, rotation=rotation, ha='center', va='top')
        # color the ticks based on the pdb
        for i, color in enumerate(colors):
            plt.gca().get_xticklabels()[i].set_color(color)
        plt.title(xlabel, fontsize=fontsize)
        plt.tight_layout()
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-heatmap-target-horizontal.png', bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-heatmap-target-horizontal.svg', format="svg", bbox_inches="tight", dpi=300, transparent=False)
        plt.close()

        # make colorbar with values of heatmap, but don't show heatmap!
        make_colorbar(heatmap, xlabel, f'/{args.system}/plots/{args.system_name_in_csv_file}-heatmap-colorbar-target-vertical.png', 'viridis', (1.1, 2), 'vertical')


        plt.figure(figsize=figsize)
        heatmap = predicted_values[:, np.newaxis].T
        plt.imshow(heatmap, aspect='auto', cmap='viridis', alpha=1)
        plt.yticks([])
        plt.xticks(range(len(mutants)), mutants, fontsize=fontsize, rotation=rotation, ha='center', va='top')
        # color the ticks based on the pdb
        for i, color in enumerate(colors):
            plt.gca().get_xticklabels()[i].set_color(color)
        plt.title(ylabel, fontsize=fontsize)
        plt.tight_layout()
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-heatmap-horizontal.png', bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-heatmap-horizontal.svg', format="svg", bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-heatmap-horizontal.pdf', format="pdf", bbox_inches="tight", dpi=300, transparent=False)
        plt.close()

        # make colorbar with values of heatmap, but don't show heatmap!
        make_colorbar(heatmap, ylabel, f'/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-heatmap-colorbar-vertical.png', 'viridis', (0.9, 2), 'vertical')


    else:
        # make 2D single-point-mutation heatmap!

        wt_seq = system_name_to_wt_seq(args.system_name_in_csv_file)

        sequences = df['sequence'].values
        targets = df[target_column].values
        predictions = df[prediction_column].values

        mutants, num_mutants = zip(*[compute_mutant(wt_seq, seq) for seq in sequences])
        mutants = np.array(mutants)
        num_mutants = np.array(num_mutants)

        assert np.all(np.logical_or(num_mutants == 1, num_mutants == 0))

        heatmap_target = np.full((len(wt_seq), 20), np.nan)
        heatmap_predicted = np.full((len(wt_seq), 20), np.nan)

        wt_coords = []

        for i, (mutant, target, prediction) in enumerate(zip(mutants, targets, predictions)):

            if len(mutant) == 0: # WT!
                for j, aa in enumerate(wt_seq):
                    heatmap_target[j, AA_TO_IND[aa]] = target
                    heatmap_predicted[j, AA_TO_IND[aa]] = prediction
                    wt_coords.append((j, AA_TO_IND[aa]))
            else:
                aa_wt = mutant[0]
                aa_mt = mutant[-1]
                pos = int(mutant[1:-1]) - 1
                heatmap_target[pos, AA_TO_IND[aa_mt]] = target
                heatmap_predicted[pos, AA_TO_IND[aa_mt]] = prediction
        
        # make the heatmaps
        figsize = (4, 6)
        fontsize = 20

        heatmap_target = heatmap_target.T
        heatmap_predicted = heatmap_predicted.T

        plt.figure(figsize=figsize)
# 
        limits = (-9, 8)
        plt.imshow(heatmap_target, aspect='auto', cmap='viridis', alpha=1)#, vmin=limits[0], vmax=limits[1])
        # plot boxes around the wildtype residues
        for wt_coord in wt_coords:
            plt.plot([wt_coord[0]-0.5, wt_coord[0]+0.5, wt_coord[0]+0.5, wt_coord[0]-0.5, wt_coord[0]-0.5], [wt_coord[1]-0.5, wt_coord[1]-0.5, wt_coord[1]+0.5, wt_coord[1]+0.5, wt_coord[1]-0.5], c='red', lw=2.0)
        plt.xticks(range(len(wt_seq)), wt_seq, fontsize=fontsize)
        plt.yticks(range(20), IND_TO_AA, fontsize=fontsize)
        plt.tight_layout()
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-heatmap-target.png', bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-heatmap-target.svg', format="svg", bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-heatmap-target.pdf', format="pdf", bbox_inches="tight", dpi=300, transparent=False)
        plt.close()
        make_colorbar(heatmap_target, xlabel, f'/{args.system}/plots/{args.system_name_in_csv_file}-heatmap-colorbar-target.png', 'viridis', (3.6, 0.8)) #, limits=limits)
        make_colorbar(heatmap_target, xlabel, f'/{args.system}/plots/{args.system_name_in_csv_file}-heatmap-colorbar-target-vertical.png', 'viridis', (1.2, 3.0), orientation='vertical') #, limits=limits)

        plt.figure(figsize=figsize)
        # limits = (12, 56)
        plt.imshow(heatmap_predicted, aspect='auto', cmap='viridis', alpha=1)#, vmin=limits[0], vmax=limits[1])
        # plot boxes around the wildtype residues
        for wt_coord in wt_coords:
            plt.plot([wt_coord[0]-0.5, wt_coord[0]+0.5, wt_coord[0]+0.5, wt_coord[0]-0.5, wt_coord[0]-0.5], [wt_coord[1]-0.5, wt_coord[1]-0.5, wt_coord[1]+0.5, wt_coord[1]+0.5, wt_coord[1]-0.5], c='red', lw=2.0)
        plt.xticks(range(len(wt_seq)), wt_seq, fontsize=fontsize)
        plt.yticks(range(20), IND_TO_AA, fontsize=fontsize)
        plt.tight_layout()
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-heatmap.png', bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-heatmap.svg', format="svg", bbox_inches="tight", dpi=300, transparent=False)
        plt.savefig(f'/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-heatmap.pdf', format="pdf", bbox_inches="tight", dpi=300, transparent=False)
        plt.close()
        make_colorbar(heatmap_predicted, ylabel, f'/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-heatmap-colorbar.png', 'viridis', (3.6, 0.8)) #, limits=limits)
        make_colorbar(heatmap_predicted, ylabel, f'/{args.system}/plots/{args.system_name_in_csv_file}-{args.model_version}-{args.with_relaxation}-heatmap-colorbar-vertical.png', 'viridis', (1.0, 3.0), orientation='vertical') #, limits=limits)


