
import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
import argparse

from scipy.special import softmax, log_softmax, logsumexp

from hermes.inference.inference_hermes import predict_from_pdbfile, load_hermes_models
from hermes.utils.protein_naming import ind_to_ol_size, ol_to_ind_size

sys.path.append('../../mutation_effects/src')
from global_constants import LOGOMAKER_COLORS

AMINO_ACIDS = list(ol_to_ind_size.keys())

def entropy(pwm, epsilon=1e-10):
    """
    Calculate the position-wise entropy of a position weight matrix (PWM).
    Robust to zeros
    """
    # Normalize the PWM
    pwm = pwm / (pwm.sum(axis=0) + epsilon)

    # Calculate the entropy
    entropy = -np.sum(pwm * np.log((pwm + epsilon)), axis=0)
    return entropy


def get_pwm(args, models, hparams, finetuning_hparams):

    region_ids = [(args.chain, resnum, ' ') for resnum in list(range(args.peptide_resnum_start, args.peptide_resnum_start + args.peptide_length))]
    requested_regions = {'peptide': region_ids}

    # compute pnE and pnlogp on the structure, save in new column of csv file
    ensemble_predictions_dict = predict_from_pdbfile(args.pdbpath, models, hparams, 64, finetuning_hparams=finetuning_hparams, regions=requested_regions)

    ensemble_predictions_dict = ensemble_predictions_dict['peptide']
    
    if args.ensemble_at_logits_level:
        ensembled_peptide_pes = np.mean(ensemble_predictions_dict['logits'], axis=0)
        ensembled_peptide_probabilities = softmax(ensembled_peptide_pes, axis=1)
    else:
        ensembled_peptide_probabilities = np.mean(ensemble_predictions_dict['probabilities'], axis=0)

    return ensembled_peptide_probabilities


def make_logoplot(pwm, output_path, use_entropy=False):

    fontsize = 18

    fig, ax = plt.subplots(1, 1, figsize=(5.4, 1.8))

    if use_entropy:
        ics = []
        for row in pwm:
            ics.append(np.log2(20) - entropy(row))
        ics = np.array(ics)
        pwm = pwm * ics[:, np.newaxis]

    # plot pwm as logoplot
    df_pwm = pd.DataFrame(pwm, columns=AMINO_ACIDS)
    logomaker.Logo(df_pwm, ax=ax, color_scheme=LOGOMAKER_COLORS)

    if use_entropy:
        ax.set_ylim([0, np.log2(20)])
    else:
        ax.set_ylim([0, 1])

    # ax.set_title(title + f'\n{len(seqs)} peptides', fontsize=16)

    ax.set_xticks(np.arange(pwm.shape[0]))
    ax.set_xticklabels(np.arange(pwm.shape[0])+1)

    ax.tick_params(axis='x', labelsize=16)

    # Remove x-tick marks but keep the labels
    ax.tick_params(axis='x', which='both', length=0)

    # ax.set_xticks([])

    # Remove all the spines
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    if not use_entropy:
        ax.spines['right'].set_visible(False)
        ax.set_yticks([])
    else:
        ax.yaxis.set_ticks_position('right')  # Set ticks on the right
        ax.yaxis.set_label_position('right')  # Set the label on the right
        ax.set_ylabel('bits', fontsize=fontsize)
        ax.set_yticks([1, 3])
        ax.tick_params(axis='y', labelsize=fontsize-2)

    if use_entropy:
        use_entropy_str = '_norm_entropy'
    else:
        use_entropy_str = ''

    plt.tight_layout()
    plt.savefig(output_path.replace('.png', f'{use_entropy_str}.png'), bbox_inches='tight')
    plt.show()



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--hermes_path', type=str, default='/gscratch/spe/gvisan01/hermes/')
    parser.add_argument('--model_version', type=str, required=True)

    parser.add_argument('--pdbpath', type=str, required=True)
    parser.add_argument('--chain', type=str, required=True)
    parser.add_argument('--peptide_length', type=int, required=True)

    parser.add_argument('--output_path_no_extension', type=str, required=True)

    parser.add_argument('--peptide_resnum_start', type=int, default=1,
                            help='Sometimes, this is not one. This might happen for example with certain structures generated in-silico.')

    parser.add_argument('--ensemble_at_logits_level', type=int, default=1, choices=[0, 1])
    
    
    args = parser.parse_args()

    # load HERMES models
    trained_models_path = os.path.join(args.hermes_path, 'trained_models', args.model_version)
    model_dir_list = [os.path.join(trained_models_path, model_rel_path) for model_rel_path in os.listdir(trained_models_path)]
    models, hparams, finetuning_hparams = load_hermes_models(model_dir_list)

    # get pwm
    pwm = get_pwm(args, models, hparams, finetuning_hparams)

    # save pwm
    df_pwm = pd.DataFrame(pwm, columns=AMINO_ACIDS)
    df_pwm.to_csv(args.output_path_no_extension + '.csv', index=True)

    # save pwm as png
    make_logoplot(pwm, args.output_path_no_extension + '.png', use_entropy=True)

