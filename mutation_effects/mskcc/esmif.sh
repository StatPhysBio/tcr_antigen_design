
for tcr in '1 2 3 4 5 6 7'
    do

    python -u ../src/score_peptides_with_esmif.py \
                --pdbdir ./pdbs \
                --output_csv_filepath './esmif/mskcc_tcr'$tcr'_ec50_sat_mut_af3-esmif-use_mt_structure=0.csv' \
                --csv_filepath 'mskcc_tcr'$tcr'_ec50_sat_mut_af3.csv' \
                --pdb_column pdb \
                --chain_column chain \
                --peptide_column sequence
    
    
done
