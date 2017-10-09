import pandas.io.sql as psql
import src.utils.python.util as _utils
import numpy as np
import pandas as pd
import logging
import sys

logger = logging.getLogger(__name__)


def process_features(df):
    """Processes mutation consequence types from probabilistic 20/20.
    """
    # rename column headers
    rename_dict = {'silent snv': 'silent'}
    df = df.rename(columns=rename_dict)

    # get nonsilent/silent
    nonsilent_to_silent = (df['non-silent snv']+df['inframe indel']+df['frameshift indel']).astype(float)/(df['silent']+1)

    ###
    # process score information
    ###
    # process conservation info for missense mutations
    if 'Total Missense MGAEntropy' in df.columns:
        num_mis = df['missense'].copy()
        #num_mis[num_mis==0] = 1
        df['Mean Missense MGAEntropy'] = np.nan
        df.loc[num_mis!=0, 'Mean Missense MGAEntropy'] = df['Total Missense MGAEntropy'][num_mis!=0] / num_mis[num_mis!=0]
        df['Mean Missense MGAEntropy'] = df['Mean Missense MGAEntropy'].fillna(df['Mean Missense MGAEntropy'].max())
        del df['Total Missense MGAEntropy']
    # calculate the mean VEST score
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

    return df


def label_gene(gene,
               oncogene=True,
               tsg=True,
               kind='onco_tsg'):
    """Label a gene according to list of oncogenes
    and tsg."""
    # set integer representation of classes
    other_num = _utils.other_label
    if oncogene: onco_num = _utils.onco_label
    if tsg: tsg_num = _utils.tsg_label if oncogene else _utils.onco_label
    smg_num = 1

    # classify genes
    if kind == 'onco_tsg':
        if gene in _utils.oncogene_set:
            return onco_num
        elif gene in _utils.tsg_set:
            return tsg_num
        else:
            return other_num
    elif kind == 'smg':
        if gene in _utils.smg_list:
            return smg_num
        else:
            return other_num


def randomize(df, prng=None):
    """Randomly shuffles the features and labels the "true" classes.

    Calls random_sort to do the random shuffling.

    Parameters
    ----------
    df : pd.DataFrame
        contains features for training
    prng : np.random.RandomState
        pseudo random number generator state

    Returns
    -------
    x : pd.DataFrame
        training features
    y : pd.Series
        true class labels
    """
    x = random_sort(df, prng)  # randomly sort data
    y = x.index.to_series().apply(label_gene)  # get gene labels
    return x, y


def random_sort(df, prng=None):
    """Randomly shuffle a DataFrame.

    NOTE: if the training data is not randomly shuffled, then
    supervised learning may find artifacts related to the order
    of the data.

    Parameters
    ----------
    df : pd.DataFrame
        dataframe with feature information

    Returns
    -------
    df : pd.DataFrame
        Randomly shuffled data frame
    """
    # get new random state if not specified
    if prng is None:
        prng = np.random.RandomState()

    # get random order
    random_indices = prng.choice(df.index.values,  # sample from 'genes'
                                 len(df),  # number of samples
                                 replace=False)  # sample without replacement

    # change order of df
    random_df = df.ix[random_indices].copy()

    return random_df


def check_num_classes(class_labels):
    class_cts = class_labels.value_counts()
    num_sufficient_uniq = (class_cts > 5).sum()
    if num_sufficient_uniq < 3:
        sys.exit('ERROR: There were either no or very few mutated oncogenes or tumor suppressor genes '
                 'found in your data! Did you supply a full pan-cancer dataset? '
                 'Or have you modified the training list of oncogenes or '
                 'tumor suppressor genes? Or did you subset your mutations to not include oncogenes/tumor suppressor genes in the training list?')


############################################
# Old feature processing functions
############################################
#import src.data_analysis.python.feature_matrix as fmat
#import src.data_analysis.python.position_entropy as pentropy

def retrieve_gene_features(conn, opts,
                           get_entropy=True):
    """Retrieve gene information from the gene_features table.

    See the gene_features module to understand the gene_features
    database table.

    Parameters
    ----------
    conn : mysql/sqlite connection
        connection to db with gene_features table
    options : dict
        options for getting info
    get_entropy : bool
        option to togle the use of entropy features.
        Since entropy features are read from a file in this function, it may
        induce a not necessary dependency on previously running commands.
        To avoid this, set get_entropy=False and then compute entropy features
        separately.

    Returns
    -------
    df : pd.dataframe
        dataframe of gene lengths
    """
    logger.info('Retrieving features of genes . . .')

    selected_cols = ['gene']

    # retrieve more features if specified by command line
    if opts['gene_length']:
        selected_cols.append('gene_length')
    if opts['mutation_rate']:
        selected_cols.append('noncoding_mutation_rate')
    if opts['replication_time']:
        selected_cols.append('replication_time')
    if opts['expression']:
        selected_cols.append('expression_CCLE as expression')
    if opts['hic']:
        selected_cols.append('HiC_compartment')
    if opts['betweeness']:
        selected_cols.append('gene_betweeness')
    if opts['degree']:
        selected_cols.append('gene_degree')

    # get info from gene_features table
    logger.info('Retrieving gene feature information from gene_features table . . . ')
    sql = "SELECT %s FROM gene_features" % ', '.join(selected_cols)
    df = psql.frame_query(sql, conn)
    df = df.set_index('gene')
    df['gene'] = df.index
    logger.info('Finished retrieving gene features from gene_features table.')

    # fill graph stats with zeros if gene not in Biogrid
    if 'gene_betweeness' in df.columns:
        df['gene_betweeness'] = df['gene_betweeness'].fillna(0)
    if 'gene_degree' in df.columns:
        df['gene_degree'] = df['gene_degree'].fillna(0)

    # get position entropy features
    if get_entropy:
        entropy_cfg = _utils.get_output_config('position_entropy')
        mutation_pos_entropy = pd.read_csv(_utils.result_dir + entropy_cfg['mutation_pos_entropy'],
                                           sep='\t', index_col=0)
        missense_pos_entropy = pd.read_csv(_utils.result_dir + entropy_cfg['missense_pos_entropy'],
                                           sep='\t', index_col=0)
        #df['mutation position entropy'] = mutation_pos_entropy['mutation position entropy']
        #df['pct of uniform mutation entropy'] = mutation_pos_entropy['pct of uniform mutation entropy']
        df['missense position entropy'] = missense_pos_entropy['missense position entropy']
        df['pct of uniform missense entropy'] = missense_pos_entropy['pct of uniform missense entropy']

    return df


def wrapper_retrieve_gene_features(opts):
    """Wrapper arround the retrieve_gene_features function in the
    features module.

    Parameters
    ----------
    opts : dict
        command line options

    Returns
    -------
    additional_features : pd.DataFrame

    """
    # get additional features
    db_cfg = _utils.get_db_config('2020plus')
    conn = sqlite3.connect(db_cfg['db'])
    additional_features = retrieve_gene_features(conn, opts, get_entropy=False)
    conn.close()
    return additional_features


def _filter_rows(df, min_ct=0):
    """Filter out rows with counts less than the minimum."""
    row_sums = df.T.sum()
    filtered_df = df[row_sums >= min_ct]
    return filtered_df


def normalize_mutational_features(df, min_count):
    """Normalizes mutation type counts and aggregate counts into
    recurrent vs deleterious.

    Parameters
    ----------
    df : pd.DataFrame
        data frame with gene names and mutation counts for each type
    min_count : int
        minimum number of mutations for a gene to be used
    """
    if 'gene' in df.columns:
        df = df.set_index('gene')  # hack to prevent dividing genes by a number
    df = _filter_rows(df, min_ct=min_count)  # drop rows below minimum total mutations
    recurrent_mutation = df['recurrent missense'] # + df['recurrent indel']
    #deleterious_mutation = df['lost stop'] + df['nonsense'] + df['frame shift'] + df['no protein'] + df['splicing mutation']
    deleterious_mutation = df['Nonstop_Mutation+Translation_Start_Site'] + df['Nonsense_Mutation'] + df['Frame_Shift_Indel'] + df['Splice_Site']
    missense_to_silent = df['Missense_Mutation'] / (df['Silent']+1).astype(float)
    row_sums = df.sum(axis=1).astype(float)
    df = df.div(row_sums-recurrent_mutation, axis=0)  # normalize each row
    df['recurrent count'] = recurrent_mutation
    df['deleterious count'] = deleterious_mutation
    #df['total'] = row_sums
    df['gene'] = df.index  # get back the gene column from the index
    return df


def generate_features(mutation_df, opts,
                      covariate_features=None):
    """Main function that generates features.

    Parameters
    ----------
    mutation_df : pd.DataFrame
        data frame containing mutations from maf/sqlitedb
    opts : dict
        dictionary containing the command line options.
        opts is not necessary if covariate_features
        is specified.
    covariate_features : pd.DataFrame (Default: None)
        if covariate data frame already obtained, then utilize
        that as input. Otherwise, retreive from sqlite database.

    Returns
    -------
    all_features : pd.DataFrame
        features with both mutational type features and covariate features
    """
    if type(covariate_features) is not pd.DataFrame:
        covariate_features = wrapper_retrieve_gene_features(opts)
    mutational_features = process_mutational_features(mutation_df)
    all_features = pd.merge(mutational_features, covariate_features,
                            how='left', on='gene')
    return all_features


def process_mutational_features(mydf):
    """Performs feature processing pipeline.

    Parameters
    ----------
    mydf : pd.DataFrame
        data frame containing the desired raw data for computation of
        features for classifier

    Returns
    -------
    proc_feat_df: pd.DataFrame
        dataframe consisting of features for classification
    """
    # rename process of columns to ensure compatability with previously
    # written code
    mydf = mydf.rename(columns={'Protein_Change': 'AminoAcid',
                                'DNA_Change': 'Nucleotide'})

    # process features
    feat_list = fmat.generate_feature_matrix(mydf, 2)
    headers = feat_list.pop(0)  # remove header row
    feat_df = pd.DataFrame(feat_list, columns=headers)  # convert to data frame
    proc_feat_df = normalize_mutational_features(feat_df, 0)
    miss_ent_df = pentropy.missense_position_entropy(mydf[['Gene', 'AminoAcid']])
    # mut_ent_df = pentropy.mutation_position_entropy(mydf[['Gene', 'AminoAcid']])

    # encorporate entropy features
    #proc_feat_df['mutation position entropy'] = mut_ent_df['mutation position entropy']
    #proc_feat_df['pct of uniform mutation entropy'] = mut_ent_df['pct of uniform mutation entropy']
    proc_feat_df['missense position entropy'] = miss_ent_df['missense position entropy']
    proc_feat_df['pct of uniform missense entropy'] = miss_ent_df['pct of uniform missense entropy']
    return proc_feat_df
