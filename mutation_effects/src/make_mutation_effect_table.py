
import os
import json
import numpy as np
import pandas as pd

MODELS = ['BLOSUM62', 'Luksza et al. $C/d$', 'Luksza et al. $C$', 'TCRdock', 'TCRdock benchmark', 'TULIP', 'TAPIR', 'NetTCR-2.2', 'ESM-IF1', 'ProteinMPNN 0.02 Pep 1-by-1 Masked', 'ProteinMPNN 0.20 Pep 1-by-1 Masked', 'ProteinMPNN 0.02 Pep Full Masked', 'ProteinMPNN 0.20 Pep Full Masked', 'ProteinMPNN 0.02 TCR Full Masked', 'ProteinMPNN 0.20 TCR Full Masked', 'ProteinMPNN 0.02 Pep and TCR Full Masked', 'ProteinMPNN 0.20 Pep and TCR Full Masked', 'HERMES-fixed 0.00', 'HERMES-fixed 0.50', 'HERMES-relaxed 0.00', 'HERMES-relaxed 0.50', 'HERMES-relaxed-min_energy 0.00', 'HERMES-relaxed-min_energy 0.50']

# DATA = [('NY-ESO', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/nyeso/results/nyeso_peptide_kd_closest_metrics.json'),
#         ('TAX', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/tax/results/tax_peptide_kd_closest_metrics.json'),
#         ('MART', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mart/results/mart_peptide_kd_closest_metrics.json'),
#         ('Hsiue', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/hsiue_et_al/results/hsiue_et_al_H2_sat_mut_metrics.json'),
#         ('Luksza TCR1', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc/results/mskcc_tcr1_ec50_sat_mut_af3_metrics.json'),
#         ('Luksza TCR2', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc/results/mskcc_tcr2_ec50_sat_mut_af3_metrics.json'),
#         ('Luksza TCR3', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc/results/mskcc_tcr3_ec50_sat_mut_af3_metrics.json'),
#         ('Luksza TCR4', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc/results/mskcc_tcr4_ec50_sat_mut_af3_metrics.json'),
#         ('Luksza TCR5', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc/results/mskcc_tcr5_ec50_sat_mut_af3_metrics.json'),
#         ('Luksza TCR6', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc/results/mskcc_tcr6_ec50_sat_mut_af3_metrics.json'),
#         ('Luksza TCR7', '/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc/results/mskcc_tcr7_ec50_sat_mut_af3_metrics.json')]

DATA = [('NY-ESO', '../nyeso/results/nyeso_peptide_kd_closest_metrics.json'),
        ('TAX', '../tax/results/tax_peptide_kd_closest_metrics.json'),
        # ('MART', '../mart/results/mart_peptide_kd_closest_metrics.json'),
        ('Hsiue', '../hsiue_et_al/results/hsiue_et_al_H2_sat_mut_metrics.json'),
        ('Luksza TCR1', '../mskcc/results/mskcc_tcr1_ec50_sat_mut_af3_metrics.json'),
        ('Luksza TCR2', '../mskcc/results/mskcc_tcr2_ec50_sat_mut_af3_metrics.json'),
        ('Luksza TCR3', '../mskcc/results/mskcc_tcr3_ec50_sat_mut_af3_metrics.json'),
        ('Luksza TCR4', '../mskcc/results/mskcc_tcr4_ec50_sat_mut_af3_metrics.json'),
        ('Luksza TCR5', '../mskcc/results/mskcc_tcr5_ec50_sat_mut_af3_metrics.json'),
        ('Luksza TCR6', '../mskcc/results/mskcc_tcr6_ec50_sat_mut_af3_metrics.json'),
        ('Luksza TCR7', '../mskcc/results/mskcc_tcr7_ec50_sat_mut_af3_metrics.json')]

# make spearman table, with p-values in parentheses
df = pd.DataFrame(index=MODELS, columns=[d[0] for d in DATA])
for model in MODELS:
    for name, path in DATA:
        with open(path) as f:
            metrics = json.load(f)
        
        if model not in metrics['Spearman r']:
            continue

        spearman_r = metrics['Spearman r'][model][0]
        spearman_r_pval = metrics['Spearman r'][model][1]
        df.loc[model, name] = f'{spearman_r:.2f} ({spearman_r_pval:.2f})'

df.to_csv('../spearman_table.csv')


# make auroc tables
for metric_name, metric in zip(['auroc_grey_vs_colored', 'auroc_reliable', 'auroc_all', 'auroc_grey_vs_above_wt'],
                               ['AUROC grey vs colored', 'AUROC reliable', 'AUROC all', 'AUROC grey vs above wt']):
    
    df = pd.DataFrame(index=MODELS, columns=[d[0] for d in DATA])
    for model in MODELS:
        for name, path in DATA[-7:]: # only luksza data
            with open(path) as f:
                metrics = json.load(f)
            
            if model not in metrics[metric]:
                continue

            value = metrics[metric][model]
            df.loc[model, name] = f'{value:.2f}'

    df.to_csv(f'../{metric_name}_table.csv')



