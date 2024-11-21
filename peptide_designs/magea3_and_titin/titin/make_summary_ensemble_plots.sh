

cd /gscratch/spe/gvisan01/protein_holography-pytorch/runtime/tcr_pmhc_experiments/peptide_generation_with_sim_anneal

input_dir=/gscratch/spe/gvisan01/peptide_designs/magea3_and_titin/titin/hcnn_pyrosetta_annealing/annealing_runs
output_dir=/gscratch/spe/gvisan01/peptide_designs/magea3_and_titin/titin/hcnn_pyrosetta_annealing/plots

pdbs='5bs0.pdb.human.MH1.A-01.A.C.DE'
model_names='so3_convnet_base_ensemble so3_convnet_noise=0p5'

num_iters='100'
annealing_schedule='multiplicative'
energy_to_use='pnE'
region_to_optimize_energy_of='pocket' # 'pocket peptide complex'
dry_or_wet='dry'

for pdb in $pdbs
    do
    for model_name in $model_names
        do
        python ensemble.py \
                    --input_dir $input_dir \
                    --output_dir $output_dir \
                    --pdb $pdb \
                    --model_name $model_name \
                    --num_iters $num_iters \
                    --annealing_schedule $annealing_schedule \
                    --energy_to_use $energy_to_use \
                    --region_to_optimize_energy_of $region_to_optimize_energy_of \
                    --dry_or_wet $dry_or_wet
    done
done

