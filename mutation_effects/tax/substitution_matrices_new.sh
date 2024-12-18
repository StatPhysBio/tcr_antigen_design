

pdbs='1ao7 1qse 1qsf closest'

base_dir='./'
output_dir=$base_dir'results/'

for pdb in $pdbs
    do

    echo $pdb


    python -u ../../src/score_with_sub_matrix.py \
                --input_csv_filepath 'tax_peptide_kd_'$pdb'.csv' \
                --output_csv_filepath $output_dir'/blosum62/tax_peptide_kd_'$pdb'-blosum62-use_mt_structure=0.csv' \
                --wt_peptide_column sequence \
                --mt_peptide_column wt_seq

    python -u pretty_plots.py \
                --system tax \
                --system_name_in_csv_file 'tax_peptide_kd_'$pdb \
                --target_column'=-log10(Kd)' \
                --model blosum62 \
                --model_instance blosum62 \
                --show_wt_lines both_from_df

done
