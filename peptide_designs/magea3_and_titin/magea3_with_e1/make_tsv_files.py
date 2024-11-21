
import os
import pandas as pd

if __name__ == '__main__':

    tsvfiles = ['wildtype/wildtype.tsv',
                'blosum62/sample_peptides_from_blosum62__temperature=1.0.tsv',
                'blosum62/sample_peptides_from_blosum62__temperature=2.0.tsv',
                'hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_base_ensemble.tsv',
                'hcnn_fixed_structure/hcnn_peptides_from_fixed_structure_so3_convnet_noise=0p5.tsv',
                'hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_base_ensemble.tsv',
                'hcnn_pyrosetta_annealing/hcnn_plus_pyrosetta_annealing_peptides_so3_convnet_noise=0p5.tsv',
                'mhc_pwm/mhc_motif_peptides.tsv',
                'proteinmpnn_v_48_002_0.7/proteinmpnn_samples_v_48_002_0.7.tsv',
                'proteinmpnn_v_48_020_0.7/proteinmpnn_samples_v_48_020_0.7.tsv']

    for tsvfile in tsvfiles:

        # make E1 mutation in all peptides in the tsv file

        directory = os.path.dirname(tsvfile)
        os.makedirs(directory, exist_ok=True)

        e1_peptides = []
        df = pd.read_csv(os.path.join('../magea3', tsvfile), sep='\t')
        for peptide in df['peptide']:
            e1_peptide = 'E' + peptide[1:]
            e1_peptides.append(e1_peptide)
        df['peptide'] = e1_peptides
        df.to_csv(tsvfile, sep='\t', index=False)
