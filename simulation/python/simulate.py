import utils.python.util as _utils
import features.python.features as features
from classify.python.r_random_forest_clf import RRandomForest
from classify.python.random_forest_clf import RandomForest
from bootstrap import Bootstrap
import plot_data
import pandas as pd
import numpy as np
from scipy import stats
import sqlite3


def r_random_forest(df, opts):
    rrclf = RRandomForest(df,
                          total_iter=1,
                          other_sample_ratio=opts['other_ratio'],
                          driver_sample=opts['driver_rate'],
                          ntrees=opts['ntrees'],
                          min_ct=opts['min_count'])
    rrclf.kfold_validation()  # run random forest
    results = pd.DataFrame({'precision': [rrclf.onco_precision,
                                          rrclf.tsg_precision],
                            'recall': [rrclf.onco_recall,
                                       rrclf.tsg_recall]},
                           index=['oncogene', 'tsg'])
    return results


def py_random_forest(df, opts):
    rclf = RandomForest(df,
                        total_iter=1,
                        ntrees=opts['ntrees'],
                        min_ct=opts['min_count'])
    rclf.kfold_validation()  # run random forest
    results = pd.DataFrame({'precision': [rclf.onco_precision,
                                          rclf.tsg_precision],
                            'recall': [rclf.onco_recall,
                                       rclf.tsg_recall]},
                           index=['oncogene', 'tsg'])
    return results


def retrieve_gene_features(opts):
    """Wrapper arround the retrieve_gene_features function in the
    features module."""
    # get additional features
    db_cfg = _utils.get_db_config('genes')
    conn = sqlite3.connect(db_cfg['db'])
    additional_features = features.retrieve_gene_features(conn, opts)
    conn.close()
    return additional_features


def merge_feature_df(count_features,
                     additional_features):
    # merge the features into one data frame
    all_features = pd.merge(count_features, additional_features,
                            how='left', on='gene')
    all_features = all_features.set_index('gene')

    # re-order columns
    #cols = all_features.columns.tolist()
    #new_order = ['gene'] + cols[:cols.index('gene')] + cols[cols.index('gene')+1:]
    #all_features = all_features[new_order]  # make the gene name the first column

    return all_features


def calculate_sem(wp):
    """Calculates the standard error of the mean for a pd.Panel object.

    **Note:** The pd.Panel.apply method seems to have a bug preventing
    me from using it. So instead I am using the numpy apply function
    for calculating sem.

    **Parameters**

    wp : pd.Panel
        panel that stratifies samples

    **Returns**

    tmp_sem : pd.DataFrame
        standard error of the mean calculated along the sample axis
    """
    tmp_sem_matrix = np.apply_along_axis(stats.sem, 0, wp.values)[0]  # hack because pandas apply method has a bug
    tmp_sem = pd.DataFrame(tmp_sem_matrix,
                           columns=['precision', 'recall'],
                           index=['oncogene', 'tsg'])
    return tmp_sem


def calculate_stats(result_dict):
    """Computes mean and sem of classification performance metrics.

    **Parameters**

    result_dict : dict
        dictionary with the i'th sample as the key and data frames
        with "oncogene"/"tsg" (row) classification performance metrics
        (columns) as values

    **Returns**

    result_df : pd.DataFrame
        Data frame with mean and sem of classification performance
        metrics. (rows: "oncogene"/"tsg", columns: summarized metrics)
    """
    wp = pd.Panel(result_dict)
    tmp_means = wp.mean(axis=0)
    tmp_sem = calculate_sem(wp)
    result_df = pd.merge(tmp_means, tmp_sem,
                         left_index=True, right_index=True,
                         suffixes=(' mean', ' sem'))
    return result_df


def main(cli_opts):
    out_opts = _utils.get_output_config('feature_matrix')

    gene_df = retrieve_gene_features(cli_opts)  # features like gene length, etc

    df = pd.read_csv(_utils.result_dir + out_opts['gene_feature_matrix'],
                     sep='\t', index_col=0)

    # iterate through each sampling rate
    r_result, py_result = {}, {}
    for sample_rate in np.linspace(.1, 3, 7):
        # bootstrap object for sub-sampling
        bs = Bootstrap(df.copy(),
                       subsample=sample_rate,
                       num_samples=cli_opts['samples'])

        # iterate through each sampled data set
        r_sim_results, py_sim_results = {}, {}  # dict to store results
        for i, bs_df in enumerate(bs.dataframe_generator()):
            # create feature matrix for bootstrap sample
            bs_df = features.process_features(bs_df, 0)
            all_df = merge_feature_df(bs_df, gene_df)  # all feature info

            # run classifiers on bootstrap sampled counts
            r_sim_results[i] = r_random_forest(all_df, cli_opts)
            py_sim_results[i] = py_random_forest(all_df, cli_opts)

        # record result for a specific sample rate
        tmp_results = calculate_stats(r_sim_results)
        r_result[sample_rate] = tmp_results
        tmp_results = calculate_stats(py_sim_results)
        py_result[sample_rate] = tmp_results

    # aggregate results for plotting
    results = {'sub-sampled random forest': pd.Panel(r_result),
               'random forest': pd.Panel(py_result)}

    # plot results of simulations
    plot_data.oncogene_precision_errorbar(results,
                                          save_path='tmp.precision.png',
                                          title='Oncogene Precision while varying DB size',
                                          xlabel='Sample rate',
                                          ylabel='Precision')
    plot_data.oncogene_recall_errorbar(results,
                                       save_path='tmp.recall.png',
                                       title='Oncogene Recall while varying DB size',
                                       xlabel='Sample rate',
                                       ylabel='Recall')
