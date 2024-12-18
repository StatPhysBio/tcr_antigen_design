

hermes_dir='/gscratch/spe/gvisan01/hermes/'

pdbs='2bnr 2bnq'

model_version_list='hermes_py_000 hermes_py_050 hermes_py_000_ft_skempi_no_tcrpmhc_ddg_bi hermes_py_050_ft_skempi_no_tcrpmhc_ddg_bi'

use_mt_structure='0'

base_dir='./'
output_dir=$base_dir'results/'

# for pdb in $pdbs
#     do

#     echo $pdb

#     for model_version in $model_version_list
#         do

#         echo $model_version

#         python -u ../../src/get_score_of_requested_peptides_fixed_structure.py \
#                     --model_version $model_version \
#                     --pdbdir ./pdbs \
#                     --output_csv_filepath $output_dir$model_version'/nyeso_peptide_kd_'$pdb'-'$model_version'-use_mt_structure=0.csv' \
#                     --csv_filepath 'nyeso_peptide_kd_'$pdb'.csv' \
#                     --pdb_column pdb \
#                     --chain_column chain \
#                     --peptide_column sequence \
#                     --peptide_resnum_start 1

#         python -u pretty_plots.py \
#                     --system nyeso \
#                     --system_name_in_csv_file 'nyeso_peptide_kd_'$pdb \
#                     --target_column'=-log10(Kd)' \
#                     --prediction_column pnE \
#                     --model hermes \
#                     --model_instance $model_version \
#                     --use_mt_structure $use_mt_structure \
#                     --show_wt_lines both_from_df
    
#     done

# done


for model_version in $model_version_list
    do

    # python -u ../../src/get_score_of_requested_peptides_fixed_structure.py \
    #             --model_version $model_version \
    #             --pdbdir ./pdbs \
    #             --output_csv_filepath $output_dir$model_version'/nyeso_peptide_kd_closest-'$model_version'-use_mt_structure=0.csv' \
    #             --csv_filepath nyeso_peptide_kd_closest.csv \
    #             --pdb_column pdb \
    #             --chain_column chain \
    #             --peptide_column sequence \
    #             --peptide_resnum_start 1

    python -u pretty_plots.py \
                --system nyeso \
                --system_name_in_csv_file nyeso_peptide_kd_closest \
                --target_column'=-log10(Kd)' \
                --model hermes \
                --model_instance $model_version \
                --prediction_column pnE \
                --use_mt_structure $use_mt_structure \
                --show_wt_lines both_from_df

    python -u pretty_plots.py \
                --system nyeso \
                --system_name_in_csv_file nyeso_peptide_kd_closest \
                --target_column'=-log10(Kd)' \
                --model hermes \
                --model_instance $model_version \
                --prediction_column pnlogp \
                --use_mt_structure $use_mt_structure \
                --show_wt_lines both_from_df

done
