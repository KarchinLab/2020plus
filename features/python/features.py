import utils.python.util as _utils
import numpy as np
import sqlite3
import pandas as pd
import pandas.io.sql as psql
import plot_data
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
    entropy_cfg = _utils.get_output_config('position_entropy')
    mutation_pos_entropy = pd.read_csv(_utils.result_dir + entropy_cfg['mutation_pos_entropy'],
                                       sep='\t', index_col=0)
    missense_pos_entropy = pd.read_csv(_utils.result_dir + entropy_cfg['missense_pos_entropy'],
                                       sep='\t', index_col=0)
    df['mutation position entropy'] = mutation_pos_entropy['mutation position entropy']
    df['pct of uniform mutation entropy'] = mutation_pos_entropy['pct of uniform mutation entropy']
    df['missense position entropy'] = missense_pos_entropy['missense position entropy']
    df['pct of uniform missense entropy'] = missense_pos_entropy['pct of uniform missense entropy']

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


def _filter_rows(df, min_ct=0):
    """Filter out rows with counts less than the minimum."""
    row_sums = df.T.sum()
    filtered_df = df[row_sums >= min_ct]
    return filtered_df


def process_features(df, min_count):
    """Process mutation type counts.

    **Parameters**

    df : pd.DataFrame
        data frame with gene names and mutation counts for each type
    min_count : int
        minimum number of mutations for a gene to be used
    """
    if 'gene' in df.columns:
        df = df.set_index('gene')  # hack to prevent dividing genes by a number
    df = _filter_rows(df, min_ct=min_count)  # drop rows below minimum total mutations
    recurrent_mutation = df['recurrent missense'] + df['recurrent indel']
    deleterious_mutation = df['lost stop'] + df['nonsense'] + df['frame shift'] + df['no protein'] + df['splicing mutation']
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
    in_opts = _utils.get_input_config('classifier')
    out_opts = _utils.get_output_config('features')
    count_opts = _utils.get_output_config('feature_matrix')
    result_opts = _utils.get_input_config('result')
    db_cfg = _utils.get_db_config('genes')

    # read in mutation counts generated in data_analysis folder
    logger.info('Getting count features . . .')
    count_features = pd.read_csv(_utils.result_dir + count_opts['gene_feature_matrix'],
                                 sep='\t')
    count_features = process_features(count_features,
                                      min_count=options['min_count'])
    logger.info('Finished getting count features.')

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

    # plot mutation histogram for olfactory receptor genes
    plot_data.or_gene_hist(all_features,
                           result_opts['feature_plot_dir'] + out_opts['or_hist'])


    # scatter plot of correlation between recurrent counts
    # and total mutation counts
    tmpdf = all_features[(all_features['recurrent count']<300) & (all_features['total']<300)]
    plot_data.correlation_plot(tmpdf,
                               'recurrent count', 'total',
                               save_path=result_opts['feature_plot_dir'] + out_opts['recur_vs_total_cor'],
                               xlabel='Recurrent Mutations',
                               ylabel='Total Mutations',
                               title='Recurrent ($<300$) vs Total ($<300$)mutations')


    # scatter plot of correlation between deleterious counts
    # and total mutation counts
    tmpdf = all_features[(all_features['deleterious count']<300) & (all_features['total']<300)]
    plot_data.correlation_plot(tmpdf,
                               'deleterious count', 'total',
                               save_path=result_opts['feature_plot_dir'] + out_opts['del_vs_total_cor'],
                               xlabel='Deleterious Mutations',
                               ylabel='Total Mutations',
                               title='Deleterious ($<300$) vs Total ($<300$)mutations')
