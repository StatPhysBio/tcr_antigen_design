

hermes_dir='/gscratch/spe/gvisan01/hermes/'

system_names='nyeso_peptide_kd_closest_af3_yes_template nyeso_peptide_kd_closest_af3_no_template'

model_version_list='hermes_py_000 hermes_py_050'
use_mt_structure='0'

base_dir='./'
output_dir=$base_dir'results/'

for system_name in $system_names
    do

    for model_version in $model_version_list
        do

        python -u ../src/get_score_of_requested_peptides_fixed_structure.py \
                    --model_version $model_version \
                    --pdbdir ./pdbs \
                    --output_csv_filepath $output_dir$model_version'/'$system_name'-'$model_version'-use_mt_structure=0.csv' \
                    --csv_filepath $system_name'.csv' \
                    --pdb_column pdb \
                    --chain_column chain \
                    --peptide_column sequence \
                    --peptide_resnum_start 1

        python -u pretty_plots.py \
                    --system nyeso \
                    --system_name_in_csv_file $system_name \
                    --target_column'=-log10(Kd)' \
                    --model hermes \
                    --model_instance $model_version \
                    --prediction_column pnE \
                    --use_mt_structure $use_mt_structure \
                    --show_wt_lines both_from_df
    
    done

done
