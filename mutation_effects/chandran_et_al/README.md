
Experiments from the paper: **Immunogenicity and therapeutic targeting of a public neoantigen derived from mutated PIK3CA**

Alanine and glycine scans for two TCRs bound to a mutated PIK3CA antigen

Originally requested by Bingxu

I hand-wrote the csv file from Figure 4b and Extended Data Figure 3d of the paper. I could not find the raw data, but this should be good enough.
- EDIT: emailed authors, got raw data

## What we should correlate with
The figures provide cytokine production (TFN\alpha) for each mutation as a *percent of the wildtype* (the figure says "% max" but I think they mean wildtype because in one instance the value is higher than 100%).
Under this paradigm, the value associated with the standard `(F_mt - F_wt)` (`F` indicates function) which we usually correlate with, is calculated as `(p - 1)F_wt` where `p` is the percent of wildtype.
For the purposes of pearson and spearman correlation, we can just correlate with `p` directly.


**PROBLEM:** the PDB files, when I open them in pymol, have the TCR and the pMHC misaligned! However, when visualizing the structures on the PDB website, they look fine. I do not understand how this could be possible.
I tried downloading the mmCIF file, same issue. I tried *exporting the model I see on the website into mmCIF*, but then I have trouble converting it to PDB. One conversion tool threw an error. Another one messes up the MHC completely.
- EDIT: asked Bingxu and figured it out. Can expand the crystal in pymol.
