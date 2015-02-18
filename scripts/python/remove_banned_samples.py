"""This script removes all samples from the mutation input
from a list provided as input.
"""

import pandas as pd
import argparse

def parse_arguments():
    info = 'Remove a list of samples from mutations'
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-b', '--ban',
                        type=str, required=True,
                        help='File with a single column with sample IDs to '
                        'remove from the mutation file (# comment lines allowed)')
    parser.add_argument('-m', '--mutations',
                        type=str, required=True,
                        help='file containing mutations as input')
    #parser.add_argument('-c', '--column',
                        #type=str, required=True,
                        #help='Either column name or column number (0-indexed)'
                        #' that contains the sample IDs to be filtered')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='output filename after filtering')

    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in sample ids to filter from data
    with open(opts['ban']) as handle:
        ban_list = [line.strip() for line in handle if not line.startswith('#')]

    # get column with sample IDs
    #try:
        #int(opts['column'])  # will fail if not int
        #is_int = True
        #col = int(opts['column'])
    #except:
        #is_int = False
        #col = opts['column']

    # read in mutations
    mut_df = pd.read_csv(opts['mutations'], sep='\t')

    # fix any potential problem with extra space for endometrial carcinoma
    mut_df['Tumor_Type'][mut_df['Tumor_Type']=='Endometrial Carcinoma '] = 'Endometrial Carcinoma'

    # remove mutations from banned sample IDs
    mut_df = mut_df[~mut_df['Tumor_Sample'].isin(ban_list)]

    # save results
    mut_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == '__main__':
    opts = parse_arguments()
    main(opts)
