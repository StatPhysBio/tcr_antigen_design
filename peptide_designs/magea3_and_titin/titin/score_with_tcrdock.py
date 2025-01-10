


import os


if __name__ == '__main__':

    curr_dir = os.path.dirname(os.path.realpath(__file__))

    # os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
    #                         -i {curr_dir}/wildtype/wildtype.tsv \
    #                         -o {curr_dir}/wildtype/tcrdock_output \
    #                         -o2 {curr_dir}/wildtype")

    # os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
    #                         -i {curr_dir}/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble.tsv \
    #                         -o {curr_dir}/hcnn_fixed_structure/tcrdock_output_base \
    #                         -o2 {curr_dir}/hcnn_fixed_structure")
    
    # os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
    #                         -i {curr_dir}/hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5.tsv \
    #                         -o {curr_dir}/hcnn_fixed_structure/tcrdock_output_noise=0p5 \
    #                         -o2 {curr_dir}/hcnn_fixed_structure")
    
    # os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
    #                         -i {curr_dir}/mhc_pwm/mhc_motif_peptides.tsv \
    #                         -o {curr_dir}/mhc_pwm/tcrdock_output \
    #                         -o2 {curr_dir}/mhc_pwm")
    
    # os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
    #                         -i {curr_dir}/proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7.tsv \
    #                         -o {curr_dir}/proteinmpnn_v_48_002_0.7/tcrdock_output \
    #                         -o2 {curr_dir}/proteinmpnn_v_48_002_0.7")

    # os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
    #                         -i {curr_dir}/proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7.tsv \
    #                         -o {curr_dir}/proteinmpnn_v_48_020_0.7/tcrdock_output \
    #                         -o2 {curr_dir}/proteinmpnn_v_48_020_0.7")

    os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
                            -i {curr_dir}/blosum62/sample_peptides_from_blosum62__temperature=1.0.tsv \
                            -o {curr_dir}/blosum62/tcrdock_output_T=1 \
                            -o2 {curr_dir}/blosum62")

    os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
                            -i {curr_dir}/blosum62/sample_peptides_from_blosum62__temperature=2.0.tsv \
                            -o {curr_dir}/blosum62/tcrdock_output_T=2 \
                            -o2 {curr_dir}/blosum62")

    os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
                            -i {curr_dir}/blosum62/sample_peptides_from_blosum62__temperature=3.0.tsv \
                            -o {curr_dir}/blosum62/tcrdock_output_T=3 \
                            -o2 {curr_dir}/blosum62")

    # os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
    #                         -i {curr_dir}/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble.tsv \
    #                         -o {curr_dir}/hcnn_pyrosetta_annealing/tcrdock_output_base \
    #                         -o2 {curr_dir}/hcnn_pyrosetta_annealing")

    # os.system(f"python -u /gscratch/spe/gvisan01/TCRdock-copy/tcrdock_pipeline_on_tsvfile.py \
    #                         -i {curr_dir}/hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5.tsv \
    #                         -o {curr_dir}/hcnn_pyrosetta_annealing/tcrdock_output_noise=0p5 \
    #                         -o2 {curr_dir}/hcnn_pyrosetta_annealing")
