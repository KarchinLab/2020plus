import argparse
import pandas as pd
import IPython


def parse_arguments():
    info = 'Format prob2020 output into summarized features'
    parser = argparse.ArgumentParser(description=info)
    help_str = 'simmulate_summary output'
    parser.add_argument('-s', '--summary',
                        type=str, required=True,
                        help=help_str)
    help_str = 'TSG output from probabilistic 20/20'
    parser.add_argument('-tsg-test', '--tsg-test',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Oncogene output from probabilistic 20/20'
    parser.add_argument('-og-test', '--og-test',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Mutsigcv covariate features'
    parser.add_argument('-c', '--covariates',
                        type=str, required=True,
                        help=help_str)
    help_str = 'File with list of oncogenes'
    parser.add_argument('-og', '--og',
                        type=str, required=True,
                        help=help_str)
    help_str = 'File with list of TSGs'
    parser.add_argument('-tsg', '--tsg',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Output feature file for 20/20+'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def process_features(df):
    """Processes mutation consequence types from probabilistic 20/20.
    """
    # rename column headers
    rename_dict = {'silent snv': 'silent'}
    df = df.rename(columns=rename_dict)

    # drop id col
    df = df.drop(['ID', 'non-silent snv'], axis=1)

    # handle mutation counts
    count_cols = ['silent', 'nonsense', 'lost stop', 'lost start', 'missense',
                  'recurrent missense', 'splice site', 'inframe indel', 'frameshift indel']
    mycounts = df[count_cols]
    df['missense'] -= df['recurrent missense']

    # calculate features
    miss_to_silent = (df['missense']+df['recurrent missense']+1).astype(float)/(df['silent']+1)

    # normalize out of total mutations
    total_cts = df[count_cols].sum(axis=1)
    norm_cts = mycounts.div(total_cts.astype(float), axis=0)
    norm_cts = norm_cts.fillna(0.0)

    # combine lost stop and lost start
    lost_start_stop = norm_cts['lost stop'] + norm_cts['lost start']
    df = df.drop(['lost start', 'lost stop'], axis=1)
    norm_cts = norm_cts.drop(['lost start', 'lost stop'], axis=1)
    count_cols.pop(3)
    count_cols.pop(2)

    df[count_cols] = norm_cts
    df['lost start/stop'] = lost_start_stop
    df['missense to silent'] = miss_to_silent

    return df


def main(opts):
    # read in prob 20/20 files
    count_df = pd.read_csv(opts['summary'], sep='\t')
    tsg_test_df = pd.read_csv(opts['tsg_test'], sep='\t')
    og_test_df = pd.read_csv(opts['og_test'], sep='\t')
    covar_df = pd.read_csv(opts['covariates'], sep='\t')
    og_test_df = og_test_df.rename(columns={'gene':'Gene'})
    tsg_test_df = tsg_test_df.rename(columns={'gene':'Gene'})

    # read in list of oncogenes / tsgs
    with open(opts['og']) as handle:
        og_list = [line.strip() for line in handle]
    with open(opts['tsg']) as handle:
        tsg_list = [line.strip() for line in handle]

    # make feature matrix
    feature_df = process_features(count_df)
    tsg_test_cols = ['Gene', 'deleterious p-value']
    feature_df = pd.merge(feature_df, tsg_test_df[tsg_test_cols],
                          how='left', on='Gene')
    og_test_cols = ['Gene', 'entropy p-value']
    feature_df = pd.merge(feature_df, og_test_df[og_test_cols],
                          how='left', on='Gene')

    # add covariate feature columns
    covar_cols = ['gene',
                  #'expression_CCLE', 'replication_time',
                  'noncoding_mutation_rate',
                  #'HiC_compartment',
                  ]
    covar_df = covar_df[covar_cols].rename(columns={'gene': 'Gene'})
    feature_df = pd.merge(feature_df, covar_df,
                          how='left', on='Gene')

    # fill na values
    rename_dict = {'Gene': 'gene'}
    feature_df = feature_df.rename(columns=rename_dict)
    feature_df = feature_df.fillna(feature_df.mean())

    feature_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
