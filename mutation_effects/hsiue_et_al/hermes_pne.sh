

hermes_dir='/gscratch/spe/gvisan01/hermes/'

model_version_list='hermes_py_000 hermes_py_050' # hermes_py_000 hermes_py_050

use_mt_structure='0'

base_dir='./'
output_dir=$base_dir'results/'

for model_version in $model_version_list
    do

    echo $model_version

    python -u ../../src/get_score_of_requested_peptides_fixed_structure.py \
                --model_version $model_version \
                --pdbdir ./pdbs \
                --output_csv_filepath $output_dir$model_version'/hsiue_et_al_H2_sat_mut-'$model_version'-use_mt_structure=0.csv' \
                --csv_filepath hsiue_et_al_H2_sat_mut.csv \
                --pdb_column wt_pdb \
                --chain_column mutant_chain \
                --peptide_column sequence \
                --peptide_resnum_start 1 \
                --same_peptide_region 1

    python -u pretty_plots.py \
                --system hsiue_et_al \
                --system_name_in_csv_file 'hsiue_et_al_H2_sat_mut' \
                --target_column'=IFN_gamma (pg/ml)' \
                --model hermes \
                --model_instance $model_version \
                --prediction_column pnE \
                --use_mt_structure $use_mt_structure \
                --show_wt_lines both_from_df

done
