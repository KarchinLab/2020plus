#!/usr/bin/env python
"""
Uses the CRAVAT Variant file from SNVGet to filter out lines
in the davoli MAF file that had mappability warnings.
"""
import pandas as pd
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()
    help_str = 'Data set from Davoli et al found at Elledge lab website'
    parser.add_argument('-c', '--cravat',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Path to MAF file from davoli2maf.py script'
    parser.add_argument('-m', '--maf',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Output MAF path after filtering using CRAVAT txt file'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def main(opts):
    cravat_df = pd.read_csv(opts['cravat'], sep='\t')
    davoli_df = pd.read_csv(opts['maf'], sep='\t')

    # filter list based on mappability warning from cravat
    prev_len = len(davoli_df)
    passed_qc = cravat_df[cravat_df['Mappability Warning'].isnull()]['ID']
    davoli_df = davoli_df.ix[passed_qc]
    after_len = len(davoli_df)
    print('Before mappability filtering: {0} lines'.format(prev_len))
    print('After mappability filtering: {0} lines'.format(after_len))
    print('Line difference: {0}'.format(prev_len-after_len))

    # save filtered file
    davoli_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
