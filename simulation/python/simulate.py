import utils.python.util as _utils
import features.python.features as features
from classify.python.r_random_forest_clf import RRandomForest
from classify.python.random_forest_clf import RandomForest
from classify.python.multinomial_nb_clf import MultinomialNaiveBayes
from classify.python.vogelstein_classifier import VogelsteinClassifier
import sklearn.metrics as metrics
from bootstrap import Bootstrap
from random_sample_names import RandomSampleNames
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
                                       rrclf.tsg_recall],
                            'ROC AUC': [rrclf.onco_mean_roc_auc,
                                        rrclf.tsg_mean_roc_auc],
                            'PR AUC': [rrclf.onco_mean_pr_auc,
                                       rrclf.tsg_mean_pr_auc]},
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
                                       rclf.tsg_recall],
                            'ROC AUC': [rclf.onco_mean_roc_auc,
                                        rclf.tsg_mean_roc_auc],
                            'PR AUC': [rclf.onco_mean_pr_auc,
                                       rclf.tsg_mean_pr_auc]},
                           index=['oncogene', 'tsg'])
    return results


def naive_bayes(df, opts):
    nbclf = MultinomialNaiveBayes(df, min_ct=opts['min_count'], total_iter=1)
    nbclf.kfold_validation()  # run random forest
    results = pd.DataFrame({'precision': [nbclf.onco_precision,
                                          nbclf.tsg_precision],
                            'recall': [nbclf.onco_recall,
                                       nbclf.tsg_recall],
                            'ROC AUC': [nbclf.onco_mean_roc_auc,
                                        nbclf.tsg_mean_roc_auc],
                            'PR AUC': [nbclf.onco_mean_pr_auc,
                                       nbclf.tsg_mean_pr_auc]},
                           index=['oncogene', 'tsg'])
    return results


def vogelstein_classification(df, opts):
    x, y_true = features.randomize(df)  # don't need randomization but gets x, y
    y_true_onco = (y_true==1).astype(int)
    y_true_tsg = (y_true==2).astype(int)

    # input for 20/20 rule classifier
    ct_list = x[['recurrent count', 'deleterious count', 'total']].values

    # predict
    vclf = VogelsteinClassifier(min_count=opts['min_count'])
    y_pred = pd.Series(vclf.predict_list(ct_list))
    y_pred_onco = (y_pred=='oncogene').astype(int)
    y_pred_tsg = (y_pred=='tsg').astype(int)
    onco_prec, onco_recall, onco_fscore, onco_support = metrics.precision_recall_fscore_support(y_true_onco, y_pred_onco)
    tsg_prec, tsg_recall, tsg_fscore, tsg_support = metrics.precision_recall_fscore_support(y_true_tsg, y_pred_tsg)

    # aggregate results
    results = pd.DataFrame({'precision': [onco_prec[1],
                                          tsg_prec[1]],
                            'recall': [onco_recall[1],
                                       tsg_recall[1]]},
                           index=['oncogene', 'tsg'])

    return results


def retrieve_gene_features(opts):
    """Wrapper arround the retrieve_gene_features function in the
    features module."""
    # get additional features
    db_cfg = _utils.get_db_config('champ')
    conn = sqlite3.connect(db_cfg['db'])
    additional_features = features.retrieve_gene_features(conn, opts)
    conn.close()

    # drop mutation entropy features
    additional_features = additional_features.drop('mutation position entropy', 1)
    additional_features = additional_features.drop('missense position entropy', 1)
    additional_features = additional_features.drop('pct of uniform mutation entropy', 1)
    additional_features = additional_features.drop('pct of uniform missense entropy', 1)
    return additional_features


def merge_feature_df(count_features,
                     additional_features):
    # merge the features into one data frame
    all_features = pd.merge(count_features, additional_features,
                            how='left', on='gene')

    # re-order columns
    #cols = all_features.columns.tolist()
    #new_order = ['gene'] + cols[:cols.index('gene')] + cols[cols.index('gene')+1:]
    #all_features = all_features[new_order]  # make the gene name the first column

    return all_features


def calculate_sem(wp, columns):
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
    tmp_sem_matrix = np.apply_along_axis(stats.sem, 0, wp.values)  # hack because pandas apply method has a bug
    tmp_sem = pd.DataFrame(tmp_sem_matrix,
                           columns=columns,
                           index=['oncogene', 'tsg'])
    return tmp_sem


def calculate_stats(result_dict,
                    metrics=['precision', 'recall', 'ROC AUC', 'PR AUC']):
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
    tmp_sem = calculate_sem(wp, metrics)
    result_df = pd.merge(tmp_means, tmp_sem,
                         left_index=True, right_index=True,
                         suffixes=(' mean', ' sem'))
    return result_df


def main(cli_opts):
    out_opts = _utils.get_output_config('feature_matrix')
    sim_opts = _utils.get_output_config('simulation')
    champ_db_path = _utils.get_db_config('champ')['db']
    conn = sqlite3.connect(champ_db_path)

    gene_df = retrieve_gene_features(cli_opts)  # features like gene length, etc

    df = pd.read_csv(_utils.result_dir + out_opts['gene_feature_matrix'],
                     sep='\t', index_col=0)

    # iterate through each sampling rate
    r_result, py_result, nb_result, v_result = {}, {}, {}, {}
    for sample_rate in np.linspace(cli_opts['lower_sample_rate'],
                                   cli_opts['upper_sample_rate'],
                                   cli_opts['num_sample_rate']):
        if cli_opts['bootstrap']:
            # bootstrap object for sub-sampling
            dfg = Bootstrap(df.copy(),
                            subsample=sample_rate,
                            num_samples=cli_opts['samples'])
        else:
            # sample with out replacement of sample names
            dfg = RandomSampleNames(sub_sample=sample_rate,
                                    num_iter=cli_opts['samples'],
                                    db_conn=conn)

        # iterate through each sampled data set
        r_sim_results, py_sim_results, nb_sim_results, v_sim_results = {}, {}, {}, {}
        for i, all_df in enumerate(dfg.dataframe_generator()):
            if cli_opts['bootstrap']:
                # create feature matrix for bootstrap sample
                all_df = features.process_features(all_df, 0)

            # merge gene features
            all_df = merge_feature_df(all_df, gene_df)  # all feature info
            all_df = all_df.set_index('gene')  # get rid of gene column

            # drop entropy columns for bootstrap sampling
            #if cli_opts['bootstrap']:
                #all_df = all_df.drop('mutation position entropy', 1)
                #all_df = all_df.drop('missense position entropy', 1)
                #all_df = all_df.drop('pct of uniform mutation entropy', 1)
                #all_df = all_df.drop('pct of uniform missense entropy', 1)

            # run classifiers on bootstrap sampled counts
            r_sim_results[i] = r_random_forest(all_df, cli_opts)
            py_sim_results[i] = py_random_forest(all_df, cli_opts)
            nb_sim_results[i] = naive_bayes(all_df, cli_opts)
            v_sim_results[i] = vogelstein_classification(all_df, cli_opts)

        # record result for a specific sample rate
        tmp_results = calculate_stats(r_sim_results)
        r_result[sample_rate] = tmp_results
        tmp_results = calculate_stats(py_sim_results)
        py_result[sample_rate] = tmp_results
        tmp_results = calculate_stats(nb_sim_results)
        nb_result[sample_rate] = tmp_results
        tmp_results = calculate_stats(v_sim_results, ['precision', 'recall'])
        v_result[sample_rate] = tmp_results

    # aggregate results for plotting
    results = {'sub-sampled random forest': pd.Panel(r_result),
               'random forest': pd.Panel(py_result),
               'naive bayes': pd.Panel(nb_result)}

    # plot results of simulations
    tmp_save_path = _utils.sim_plot_dir + sim_opts['pr_plot']
    plot_data.oncogene_pr_auc_errorbar(results,
                                       save_path=tmp_save_path,
                                       title='Oncogene PR AUC vs. DB size',
                                       xlabel='Sample rate',
                                       ylabel='PR AUC')
    tmp_save_path = _utils.sim_plot_dir + sim_opts['roc_plot']
    plot_data.oncogene_roc_auc_errorbar(results,
                                        save_path=tmp_save_path,
                                        title='Oncogene ROC AUC vs. DB size',
                                        xlabel='Sample rate',
                                        ylabel='ROC AUC')

    # since the vogelstein classifier doesn't predict probabilities
    # I can't generate a PR or ROC curve. However, I can evaluate
    # metrics like precision and recall
    results['20/20 classifier'] = pd.Panel(v_result)

    tmp_save_path = _utils.sim_plot_dir + sim_opts['precision_plot']
    plot_data.oncogene_precision_errorbar(results,
                                          save_path=tmp_save_path,
                                          title='Oncogene Precision vs. DB size',
                                          xlabel='Sample rate',
                                          ylabel='Precision')
    tmp_save_path = _utils.sim_plot_dir + sim_opts['recall_plot']
    plot_data.oncogene_recall_errorbar(results,
                                       save_path=tmp_save_path,
                                       title='Oncogene Recall vs. DB size',
                                       xlabel='Sample rate',
                                       ylabel='Recall')

