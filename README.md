# T-cell receptor specificity landscape revealed through de novo peptide design

T-cells play a key role in adaptive immunity by mounting specific responses against diverse pathogens. An effective binding between  T-cell receptors (TCRs) and pathogen-derived peptides presented on Major Histocompatibility Complexes (MHCs) mediate an immune response. However, predicting these interactions remains challenging due to limited functional data on  T-cell reactivities. Here, we introduce a computational approach to predict TCR interactions with peptides presented on MHC class I alleles, and to design novel immunogenic peptides for specified TCR-MHC complexes. Our method leverages HERMES, a structure-based, physics-guided machine learning model trained on the protein universe to predict amino acid preferences based on local structural environments. Despite no direct training on TCR-pMHC data, the implicit physical reasoning in HERMES enables us to make accurate predictions of both TCR-pMHC binding affinities and T-cell activities across diverse viral epitopes and cancer neoantigens, achieving up to 72% correlation with experimental data. Leveraging our TCR recognition model, we develop a computational protocol for {de novo} design of immunogenic peptides. Through experimental validation in three TCR-MHC systems targeting viral and cancer peptides, we demonstrate that our designs---with up to five substitutions from the native sequence---activate T-cells at success rates of up to 50%. Lastly, we use our generative framework to quantify the diversity of the peptide recognition landscape for various TCR-MHC complexes, offering key insights into T-cell specificity in both humans and mice. Our approach provides a platform for immunogenic peptide and neoantigen design, opening new computational paths for T-cell vaccine development against viruses and cancer.

![Schematic](local/schematic.png)


This repository contains data and code to generate the results in the paper "T-cell receptor specificity landscape revealed through de novo peptide design" by Visani G.M. et al.\
The code is contingent upon the installation of the following three tools:
- HERMES (https://github.com/StatPhysBio/hermes)
- ProteinMPNN (fork with code to score and design peptides https://github.com/gvisani/ProteinMPNN-copy)
- TCRdock (fork with code to score peptides https://github.com/gvisani/TCRdock-copy)


## TCR-pMHC binding affinity and T-cell activity prediction

Code and benchmarking results of TCR-pMHC binding affinity prediction across different peptides can be found in `mutation_effects/nyeso` and `mutation_effects/tax`. \
Code and benchmarking results of T-cell activity prediction across different peptides can be found in `mutation_effects/mskcc` and `mutation_effects/hsiue_et_al`. \
Code for running HERMES-relaxed to score peptides can be found in `mutation_effects/src`.

## De-novo peptide design with in-vitro validation

![Peptide design](local/design_figure.png)

If you're only looking for our HERMES-made peptide designs with in-vitro T-cell activity measurements, you can find them all in `peptide_designs/All-designs.xlsx`.

Scripts to generate peptides with different algorithms, benchmarking results (via TCRdock PAE), and post-in-vitro-experiment analysis can be found in `peptide_designs/nyeso`, `peptide_designs/ebv`, `peptide_designs/magea3_and_titin`.


## Quantifying the diversity of the TCR recognition landscape

![TCR specificity](local/entropy_figure.png)

We use the entropy of HERMES-fixed's predicted Position Weight Matrix of peptides to quantify the diversity of the TCR recognition landscape. Specifically, we compute peptide entropies conditioned on TCR-MHC structures using HERMES-fixed, and compare them to the entropies of the same peptides conditioned on the MHC only, grabbed from the MHC Motif Atlas (http://mhcmotifatlas.org/home). \
Code and data can be found in `tcr_specificity`.



## Generating peptides with HERMES-*fixed*

**Install HERMES:** Clone and follow the installation steps of HERMES (https://github.com/StatPhysBio/hermes).

**Sampling script:** You can use the script `sample_peptides_with_hermes_fixed.py` to generate peptides with HERMES-*fixed*. See the docstring below:

```
usage: sample_peptides_with_hermes_fixed.py [-h] -m MODEL_VERSION -p PDB_FILE -c PEPTIDE_CHAIN -o OUTPUT_FILE [-t TEMPERATURE] [-n NUM_SAMPLES] [-b BATCH_SIZE]

Sample peptides with HERMES fixed structure

optional arguments:
  -h, --help            show this help message and exit
  -m MODEL_VERSION, --model_version MODEL_VERSION
                        Name of HERMES model you want to use. Use `hermes_py_000` or `hermes_bp_000` for the model trained without added noise to the atoms, `hermes_py_050` or `hermes_bp_050` for the model trained with 0.50 Angstrom noise to the atoms. Models with `_py_` in the name use pyrosetta to parse protein structures and were the ones used in the paper, whereas models with `_bp_` in the name use biopython.
  -p PDB_FILE, --pdb_file PDB_FILE
                        Path to the PDB file containing the TCR-pMHC structure.
  -c PEPTIDE_CHAIN, --peptide_chain PEPTIDE_CHAIN
                        Chain identifier for the peptide in the PDB file.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to the output .txt file where sampled peptide sequences will be saved (one sequence per line).
  -t TEMPERATURE, --temperature TEMPERATURE
                        Temperature for sampling from the probability distribution (default: 1.0, lower is more conservative, higher is more diverse).
  -n NUM_SAMPLES, --num_samples NUM_SAMPLES
                        Number of peptide samples to generate (note that they are *not* guaranteed to be unique). Default: 200.
  -b BATCH_SIZE, --batch_size BATCH_SIZE
                        Batch size (number of sites) for running HERMES inference. Default: 32, no need to change this except for handling memory issues.
```

**Sampling function:** You can also directly call the function `sample_peptides_with_hermes_fixed` in `sample_peptides_with_hermes_fixed.py` from your own code, with the same arguments as above (except for the output file):

```python
from sample_peptides_with_hermes_fixed import sample_peptides_with_hermes_fixed
sequences = sample_peptides_with_hermes_fixed(model_version, pdb_file, peptide_chain, temperature, num_samples, batch_size)
```

**Time:** This is very fast, taking only a handful of seconds even for a high number of sampled peptides.


## Generating peptides with HERMES-*relaxed*

**Install HERMES:** Clone and follow the installation steps of HERMES (https://github.com/StatPhysBio/hermes). You *must* include the installation of pyrosetta as well.

**Sampling script (single peptide):** Use the script `sample_single_peptide_with_hermes_relaxed.py` to generate a *single* peptide with HERMES-*relaxed*. This iterates between HERMES-*fixed* sampling and pyrosetta mutations/relaxations within a simulated annealing loop. We emphasize that this script only samples a single peptide, and should be run multiple times with different outpout files to obtain multiple peptide samples. See the docstring below:

```
usage: sample_single_peptide_with_hermes_relaxed.py [-h] --hermes_dir HERMES_DIR -m MODEL_VERSION -p PDB_FILE -o OUTPUT_FILE -c PEPTIDE_CHAIN [--initial_sequence INITIAL_SEQUENCE] --schedule_name SCHEDULE_NAME [--iters ITERS] [--energy_to_use {pnE,pnlogp}]
                                                    [--region_to_optimize_energy_of {peptide,pocket,complex}] [--dont_relax_peptide_backbone] [--save_metric_history] [--write_pdb] [--max_atoms MAX_ATOMS] [--peptide_acceptance] [--force_best_mutations]
                                                    [--add_crystal_water] [--logging LOGGING] [--pyrosetta_logging PYROSETTA_LOGGING] [--rosetta_logging ROSETTA_LOGGING]

optional arguments:
  -h, --help            show this help message and exit
  --hermes_dir HERMES_DIR
                        Path to HERMES directory containing the `trained_models` directory with the model you want to use for inference.
  -m MODEL_VERSION, --model_version MODEL_VERSION
                        Name of HERMES model you want to use. Use `hermes_py_000` or `hermes_bp_000` for the model trained without added noise to the atoms, `hermes_py_050` or `hermes_bp_050` for the model trained with 0.50 Angstrom noise to the atoms. Models with
                        `_py_` in the name use pyrosetta to parse protein structures and were the ones used in the paper, whereas models with `_bp_` in the name use biopython.
  -p PDB_FILE, --pdb_file PDB_FILE
                        Path to the PDB file containing the TCR-pMHC structure.
  -o OUTPUT_FILE, --output_file OUTPUT_FILE
                        Path to the output .txt file where sampled peptide sequences will be saved (one sequence per line).
  -c PEPTIDE_CHAIN, --peptide_chain PEPTIDE_CHAIN
                        Chain identifier for the peptide in the PDB file.
  --initial_sequence INITIAL_SEQUENCE
                        Initial peptide sequence for annealing. Empty string for random sequence
  --schedule_name SCHEDULE_NAME
                        Annealing schedule JSON file name (without the .json extension) as found in `./temperature/schedule_files`
  --iters ITERS         Number of simulated annealing iterations. Default: 100. NOTE: the schedules in this directory are designed for 100 iterations.
  --energy_to_use {pnE,pnlogp}
  --region_to_optimize_energy_of {peptide,pocket,complex}
  --dont_relax_peptide_backbone
  --save_metric_history
                        Whether to save the metric history of the simulated annealing algorithm. It will be saved in a .npy file by the same name as `output_file`.
  --write_pdb           Whether to save the relaxed PDB files with each peptide. They will be saved in a directory identified by the same name as `output_file` and suffix `__pdbs`.
  --max_atoms MAX_ATOMS
  --peptide_acceptance
  --force_best_mutations
  --add_crystal_water
  --logging LOGGING
  --pyrosetta_logging PYROSETTA_LOGGING
  --rosetta_logging ROSETTA_LOGGING
```

**Time:** This is veeeery slow due to the pyrosetta relaxations. Each round of relaxations last approximately one minute. Multiply that by 100 iterations, and you have 1.5-2 hours to sample a single peptide. We recommend running this in paralle on many CPUs. Using GPUs to speed-up the HERMES' forward pass is unnecessary as that is *not* the bottleneck of this procedure.