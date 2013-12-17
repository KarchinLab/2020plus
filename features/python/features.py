import utils.python.util as _utils
import numpy as np
import sqlite3
import pandas as pd
import pandas.io.sql as psql
import logging

logger = logging.getLogger(__name__)

def retrieve_gene_features(conn, opts):
    """Retrieve gene information from the gene_features table.

    See the gene_features module to understand the gene_features
    database table.

    **Parameters**

    conn : mysql/sqlite connection
        connection to db with gene_features table
    options : dict
        options for getting info

    **Returns**

    df : pd.dataframe
        dataframe of gene lengths
    """
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

    # get info from gene_features table
    logger.info('Retrieving gene feature information . . . ')
    sql = "SELECT %s FROM gene_features" % ', '.join(selected_cols)
    df = psql.frame_query(sql, conn)
    # df = df.set_index('gene')
    # df['gene_length'] = df['gene_length'].fillna(df['gene_length'].mean())  # replace missing with mean length
    logger.info('Finished retrieving gene features.')
    return df


def label_gene(gene,
               oncogene=True,
               tsg=True,
               kind='vogelstein'):
    """Label a gene according to Vogelstein's list of oncogenes
    and tsg."""
    # set integer representation of classes
    other_num = 0
    if oncogene: onco_num = 1
    if tsg: tsg_num = 2 if oncogene else 1
    smg_num = 1

    # classify genes
    if kind == 'vogelstein':
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


def filter_rows(df, min_ct=0):
    """Filter out rows with counts less than the minimum."""
    row_sums = df.T.sum()
    filtered_df = df[row_sums > min_ct]
    return filtered_df


def process_features(df):
    """Process mutation type counts."""
    df = df.set_index('gene')  # hack to prevent dividing genes by a number
    recurrent_mutation = df['recurrent missense'] + df['recurrent indel']
    deleterious_mutation = df['lost stop'] + df['nonsense'] + df['frame shift'] + df['no protein']
    row_sums = df.sum(axis=1).astype(float)
    df = df.div(row_sums, axis=0)  # normalize each row
    df['recurrent count'] = recurrent_mutation
    df['deleterious count'] = deleterious_mutation
    df['total'] = row_sums
    df['gene'] = df.index  # get back the gene column from the index
    return df


def randomize(df):
    """Randomly shuffles the features and labels the "true" classes.

    Calls random_sort to do the random shuffling.

    **Parameters**

    df : pd.DataFrame
        contains features for training

    **Returns**

    x : pd.DataFrame
        training features
    y : pd.Series
        true class labels
    """
    x = random_sort(df)  # randomly sort data
    y = x.index.to_series().apply(label_gene)  # get gene labels
    return x, y


def random_sort(df):
    """Randomly shuffle a DataFrame.

    NOTE: if the training data is not randomly shuffled, then
    supervised learning may find artifacts related to the order
    of the data.

    **Parameters**

    df : pd.DataFrame
        dataframe with feature information

    **Returns**

    df : pd.DataFrame
        Randomly shuffled data frame
    """
    random_indices = np.random.choice(df.index.values,  # sample from 'genes'
                                      len(df),  # number of samples
                                      replace=False)  # sample without replacement
    random_df = df.reindex(random_indices)  # change order of df
    return random_df


def main(options):
    # get configs
    #cfg_opts = _utils.get_output_config('classifier')
    in_opts = _utils.get_input_config('classifier')
    count_opts = _utils.get_output_config('feature_matrix')
    db_cfg = _utils.get_db_config('genes')

    # read in mutation counts generated in data_analysis folder
    count_features = pd.read_csv(_utils.result_dir + count_opts['gene_feature_matrix'],
                                 sep='\t')
    count_features = process_features(count_features)

    # get additional features
    conn = sqlite3.connect(db_cfg['db'])
    additional_features = retrieve_gene_features(conn, options)
    conn.close()

    # merge the features into one data frame
    all_features = pd.merge(count_features, additional_features,
                            how='left', on='gene')

    # save features to text file
    cols = all_features.columns.tolist()
    new_order = ['gene'] + cols[:cols.index('gene')] + cols[cols.index('gene')+1:]
    all_features = all_features[new_order]  # make the gene name the first column
    all_features.to_csv(in_opts['gene_feature'], sep='\t', index=False)
