import argparse
import pandas as pd
import numpy as np
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
    help_str = 'simulate_non_silent_ratio output'
    parser.add_argument('-n', '--non-silent',
                        type=str,
                        default=None,
                        help=help_str)
    help_str = 'Mutsigcv covariate features'
    parser.add_argument('-c', '--covariates',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Output feature file for 20/20+'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def process_features(df, non_silent_df=None):
    """Processes mutation consequence types from probabilistic 20/20.
    """
    # rename column headers
    rename_dict = {'silent snv': 'silent'}
    df = df.rename(columns=rename_dict)

    # calculate the mean/std for the entire cohort
    if non_silent_df is not None:
        non_silent_df = non_silent_df.rename(columns=lambda x: x.replace(' count', ''))
        non_sil_cols = ['nonsense', 'silent', 'splice site', 'lost stop',
                        'missense', 'lost start']
        total_cts = non_silent_df[non_sil_cols].sum(axis=1)
        mycounts = non_silent_df[non_sil_cols]
        nonsil_to_silent = (df['non-silent snv']+df['inframe indel']+df['frameshift indel']).astype(float)/(df['silent']+1)
        miss_to_silent = (non_silent_df['missense']).astype(float)/(non_silent_df['silent']+1)
        norm_cts = mycounts.div(total_cts.astype(float), axis=0)
        norm_cts = norm_cts.fillna(0.0)
        lost_start_stop = norm_cts['lost stop'] + norm_cts['lost start']
        norm_cts = norm_cts.drop(['lost start', 'lost stop'], axis=1)
        norm_cts['lost start and stop'] = lost_start_stop
        norm_cts['missense to silent'] = miss_to_silent
        norm_cts['non-silent to silent'] = nonsil_to_silent
        feature_means = norm_cts.mean()
        feature_stdev = norm_cts.std()

    # get nonsilent/silent
    nonsilent_to_silent = (df['non-silent snv']+df['inframe indel']+df['frameshift indel']).astype(float)/(df['silent']+1)

    # process score information
    if 'Total Missense MGAEntropy' in df.columns:
        num_mis = df['missense'].copy()
        #num_mis[num_mis==0] = 1
        df['Mean Missense MGAEntropy'] = np.nan
        df.loc[num_mis!=0, 'Mean Missense MGAEntropy'] = df['Total Missense MGAEntropy'][num_mis!=0] / num_mis[num_mis!=0]
        df['Mean Missense MGAEntropy'] = df['Mean Missense MGAEntropy'].fillna(df['Mean Missense MGAEntropy'].max())
        del df['Total Missense MGAEntropy']
    if 'Total Missense VEST Score' in df.columns:
        sum_cols = ['Total Missense VEST Score', 'lost stop',
                    'lost start', 'splice site', 'frameshift indel', 'inframe indel',
                    'nonsense']
        all_muts = ['non-silent snv', 'silent', 'inframe indel', 'frameshift indel']
        tot_vest_score = df[sum_cols].sum(axis=1).astype(float)
        num_muts = df[all_muts].sum(axis=1).astype(float)
        df['Mean VEST Score'] = tot_vest_score / num_muts
        #df['Missense VEST Score'] = df['Total Missense VEST Score'] / num_muts
        #df['VEST normalized missense position entropy'] = df['normalized missense position entropy'] * (1.-df['Missense VEST Score'])
        del df['Total Missense VEST Score']

    # drop id col
    df = df.drop(['ID', 'non-silent snv'], axis=1)

    # handle mutation counts
    count_cols = ['silent', 'nonsense', 'lost stop', 'lost start', 'missense',
                  'recurrent missense', 'splice site', 'inframe indel', 'frameshift indel']
    mycounts = df[count_cols]
    df['missense'] -= df['recurrent missense']

    # calculate features
    miss_to_silent = (df['missense']+df['recurrent missense']).astype(float)/(df['silent']+1)

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
    df['lost start and stop'] = lost_start_stop
    df['missense to silent'] = miss_to_silent
    df['non-silent to silent'] = nonsilent_to_silent

    # normalize
    if non_silent_df is not None:
        cols = ['nonsense', 'silent', 'splice site', 'lost start and stop',
                'missense', 'missense to silent', 'non-silent to silent']
        normed = (df[cols] - feature_means) #/ feature_stdev.astype(float)
        df[normed.columns] = normed

    return df


def main(opts):
    # read in prob 20/20 files
    count_df = pd.read_csv(opts['summary'], sep='\t')
    tsg_test_df = pd.read_csv(opts['tsg_test'], sep='\t')
    og_test_df = pd.read_csv(opts['og_test'], sep='\t')
    covar_df = pd.read_csv(opts['covariates'], sep='\t')
    og_test_df = og_test_df.rename(columns={'gene':'Gene'})
    tsg_test_df = tsg_test_df.rename(columns={'gene':'Gene'})

    # use cohort simulation if provided
    if opts['non_silent']:
        non_sil_df = pd.read_csv(opts['non_silent'], sep='\t')
    else:
        non_sil_df = None

    # make feature matrix
    feature_df = process_features(count_df, non_sil_df)
    tsg_test_cols = ['Gene', 'inactivating p-value']
    feature_df = pd.merge(feature_df, tsg_test_df[tsg_test_cols],
                          how='left', on='Gene')
    og_test_cols = ['Gene', 'entropy p-value', 'vest p-value', 'combined p-value']
    feature_df = pd.merge(feature_df, og_test_df[og_test_cols],
                          how='left', on='Gene')

    # add covariate feature columns
    covar_cols = ['gene',
                  'expression_CCLE',
                  'replication_time',
                  #'noncoding_mutation_rate',
                  'HiC_compartment',
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
