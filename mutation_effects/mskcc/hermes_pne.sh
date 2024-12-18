

hermes_dir='/gscratch/spe/gvisan01/hermes/'

tcrs='6 7' # 1 2 3 4 5 6 7

model_version_list='hermes_py_000 hermes_py_050'

use_mt_structure='0'

base_dir='./'
output_dir=$base_dir'results/'

for tcr in $tcrs
    do

    echo $tcr

    for model_version in $model_version_list
        do

        echo $model_version
        
        python -u ../../src/get_score_of_requested_peptides_fixed_structure.py \
                    --model_version $model_version \
                    --pdbdir ./pdbs \
                    --output_csv_filepath $output_dir$model_version'/mskcc_tcr'$tcr'_ec50_sat_mut_af3-'$model_version'-use_mt_structure=0.csv' \
                    --csv_filepath $base_dir'mskcc_tcr'$tcr'_ec50_sat_mut_af3.csv' \
                    --pdb_column wt_pdb \
                    --chain_column mutant_chain \
                    --peptide_column sequence \
                    --peptide_resnum_start 1 \
                    --same_peptide_region 1

        python -u pretty_plots.py \
                    --system mskcc \
                    --system_name_in_csv_file 'mskcc_tcr'$tcr'_ec50_sat_mut_af3' \
                    --target_column'=- delta log_ec50_M' \
                    --model hermes \
                    --model_instance $model_version \
                    --prediction_column pnE \
                    --use_mt_structure $use_mt_structure \
                    --show_wt_lines both_from_df

    
    done
done

# '=- delta log_ec50_M' \
# '=<Theta21>' \
