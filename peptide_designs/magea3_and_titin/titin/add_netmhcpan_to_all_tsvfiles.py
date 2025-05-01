
import os


if __name__ == '__main__':

    tsvfiles = ['wildtype/wildtype_w_pae_w_blosum.tsv',
                'blosum62/sample_peptides_from_blosum62__temperature=1.0_w_pae_w_blosum.tsv',
                'blosum62/sample_peptides_from_blosum62__temperature=2.0_w_pae_w_blosum.tsv',
                'blosum62/sample_peptides_from_blosum62__temperature=3.0_w_pae_w_blosum.tsv',
                'hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble_w_pae_w_blosum.tsv',
                'hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5_w_pae_w_blosum.tsv',
                'hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble_w_pae_w_blosum.tsv',
                'hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5_w_pae_w_blosum.tsv',
                'mhc_pwm/mhc_motif_peptides_w_pae_w_blosum.tsv',
                'proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7_w_pae_w_blosum.tsv',
                'proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7_w_pae_w_blosum.tsv',
                'esmif_0.7/esmif_samples_0.7_w_pae_w_blosum.tsv']

    for tsvfile in tsvfiles:
        os.system(f'python ../../src/add_netmhcpan_to_tsvfile.py \
                                                -i {tsvfile} \
                                                -m HLA-A01:01')

