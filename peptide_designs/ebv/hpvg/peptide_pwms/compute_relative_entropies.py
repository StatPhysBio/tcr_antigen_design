
import os
import numpy as np
import pandas as pd

def relative_entropy(pwm1, pwm2, epsilon=1e-10):
    """
    Calculate the relative entropy between two position weight matrices (PWMs).
    Robust to zeros
    """
    # Normalize the PWMs
    pwm1 = pwm1 / (pwm1.sum(axis=0) + epsilon)
    pwm2 = pwm2 / (pwm2.sum(axis=0) + epsilon)

    # Calculate the relative entropy
    kl_divergence = np.sum(pwm1 * np.log((pwm1 + epsilon) / (pwm2 + epsilon)))
    return kl_divergence

models = ['hermes_py_000', 'hermes_py_050']

groups = [['3mv7', 'fold_ebv_no_template_model_0', 'fold_ebv_yes_template_model_0']]

with open('relative_entropies.txt', 'w+') as f:
    for model in models:
        f.write(f'Model: {model}\n')
        for group in groups:

            exper = group[0]
            no_temp = group[1]
            yes_temp = group[2]

            pwm_exper = pd.read_csv(f'pwm_{exper}_{model}.csv', index_col=0).values
            pwm_no_temp = pd.read_csv(f'pwm_{no_temp}_{model}.csv', index_col=0).values
            pwm_yes_temp = pd.read_csv(f'pwm_{yes_temp}_{model}.csv', index_col=0).values

            # Calculate relative entropy
            rel_entropy_no_temp = relative_entropy(pwm_exper, pwm_no_temp)
            rel_entropy_yes_temp = relative_entropy(pwm_exper, pwm_yes_temp)

            # Save the results
            f.write(f'\t{exper} vs. {no_temp}: {rel_entropy_no_temp:.3f}\n')
            f.write(f'\t{exper} vs. {yes_temp}: {rel_entropy_yes_temp:.3f}\n')
        
        f.write('\n')



