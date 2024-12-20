
import os
import numpy as np
import pandas as pd
import time
import argparse

'''

1. Adds column 'is_binder_by_netmhcpan' to the input .tsv file
2. Saves the new .tsv file with the same name, unless output_tsvfile is specified

'''

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_tsvfile', type=str, required=True)
    parser.add_argument('-o', '--output_tsvfile', type=str, default=None)
    parser.add_argument('-m', '--mhc_allele', type=str, required=True, help='MHC allele for NetMHCpan, must be in format recognized by NetMHCPan')
    parser.add_argument('-p', '--peptides_column', type=str, default='peptide')
    args = parser.parse_args()

    if args.output_tsvfile is None:
        output_tsvfile = args.input_tsvfile
    else:
        output_tsvfile = args.output_tsvfile

    df = pd.read_csv(args.input_tsvfile, sep='\t')

    peptides = df[args.peptides_column].values

    min_peptide_length = min([len(peptide) for peptide in peptides])
    max_peptide_length = max([len(peptide) for peptide in peptides])
    peptide_lengths_for_netmhcpan = ','.join([str(length) for length in range(min_peptide_length, max_peptide_length+1)])
    print(f'Peptide lengths for NetMHCpan: {peptide_lengths_for_netmhcpan}')

    ## NOTE: tempfile solution isn't working for some reason, manually making tempfile with name given by time.time() so I can be pretty sure they are unique

    tempfile_in = str(time.time())
    time.sleep(0.1)
    tempfile_out = str(time.time())

    # write peptides to temp file
    with open(tempfile_in, 'w') as f:
        f.write('\n'.join(peptides))

    # run netMHCpan
    os.system(f'netMHCpan -a {args.mhc_allele} -p {tempfile_in} -l {peptide_lengths_for_netmhcpan} -xls -xlsfile {tempfile_out}')

    # read netMHCpan output
    netmhcpan_df = pd.read_csv(tempfile_out, sep='\t', skiprows=1)

    # # add is_binder_by_netmhcpan column to the input tsv file
    # df['is_binder_by_netmhc_pan'] = netmhcpan_df['NB'].values

    # save the new tsv file
    df.to_csv(output_tsvfile, sep='\t', index=False)
    
    # remove tempfiles
    os.remove(tempfile_in)
    os.remove(tempfile_out)


