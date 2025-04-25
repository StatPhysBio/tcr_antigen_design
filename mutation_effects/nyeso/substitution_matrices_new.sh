

pdbs='2bnr 2bnq closest'

sub_matrices='blosum62 luksza_cross_reactivity luksza_cross_reactivity_without_d'

base_dir='./'
output_dir=$base_dir'results/'

for sub_matrix in $sub_matrices
    do

    for pdb in $pdbs
        do

        echo $pdb


        python -u ../src/score_with_sub_matrix.py \
                    --input_csv_filepath 'nyeso_peptide_kd_'$pdb'.csv' \
                    --output_csv_filepath $output_dir'/'$sub_matrix'/nyeso_peptide_kd_'$pdb'-'$sub_matrix'-use_mt_structure=0.csv' \
                    --wt_peptide_column sequence \
                    --mt_peptide_column wt_seq \
                    --substitution_matrix $sub_matrix \

        python -u pretty_plots.py \
                    --system nyeso \
                    --system_name_in_csv_file 'nyeso_peptide_kd_'$pdb \
                    --target_column'=-log10(Kd)' \
                    --model $sub_matrix \
                    --model_instance $sub_matrix \
                    --show_wt_lines both_from_df

    done

done
