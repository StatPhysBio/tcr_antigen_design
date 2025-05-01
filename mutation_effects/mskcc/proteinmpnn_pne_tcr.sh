

proteinmpnn_dir=/gscratch/spe/gvisan01/ProteinMPNN-copy/peptide_scoring_and_design/

base_dir='/gscratch/spe/gvisan01/tcr_pmhc/mutation_effects/mskcc/'

tcrs='1 2 3 4 5 6 7' # '1 2 3 4 5 6 7'

use_mt_structure='0'

model_version_list='v_48_002 v_48_020' # v_48_030 v_48_020 v_48_002'

for tcr in $tcrs
    do

    echo $tcr

    for model_version in $model_version_list
        do

        python $proteinmpnn_dir'score_multiple_peptides_and_or_arbitrary_residues__full_masked.py' \
                    --csv_file $base_dir'mskcc_tcr'$tcr'_ec50_sat_mut_af3.csv' \
                    --folder_with_pdbs $base_dir'pdbs' \
                    --output_dir $base_dir'results/proteinmpnn_'$model_version'_full_pep_and_tcr_masked' \
                    --model_version $model_version \
                    --pdb_column wt_pdb \
                    --chain_column mutant_chain \
                    --peptides_column sequence \
                    --requested_regions all_peptide_residues interacting_TCR_residues \
                    --num_seq_per_target 10 \
                    --batch_size 10

        python -u pretty_plots.py \
                    --base_dir ../ \
                    --system mskcc \
                    --system_name_in_csv_file 'mskcc_tcr'$tcr'_ec50_sat_mut_af3' \
                    --target_column'=-log10(EC50)' \
                    --prediction_column pnlogp \
                    --model proteinmpnn \
                    --model_instance 'proteinmpnn_'$model_version'_full_pep_and_tcr_masked' \
                    --num_seq_per_target 10 \
                    --use_mt_structure $use_mt_structure \
                    --show_wt_lines both_from_df

    done
done

