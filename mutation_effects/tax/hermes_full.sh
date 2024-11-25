
# just for sanity check purposes, which passed

hermes_dir='/gscratch/spe/gvisan01/hermes/'

python -u $hermes_dir'run_hermes_on_pdbfiles.py' \
                    -m hermes_py_000 \
                    -pd ./pdbs \
                    -pn hermes_pdbs_and_chains.txt \
                    -r logprobas logits \
                    -o hermes_py_000_output.csv

