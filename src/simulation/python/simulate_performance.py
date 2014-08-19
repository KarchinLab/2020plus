import src.utils.python.util as _utils
import src.features.python.features as features
from src.classify.python.r_random_forest_clf import RRandomForest
from src.classify.python.random_forest_clf import RandomForest
from src.classify.python.multinomial_nb_clf import MultinomialNaiveBayes
from src.classify.python.vogelstein_classifier import VogelsteinClassifier
import sklearn.metrics as metrics
from bootstrap import Bootstrap
from random_sample_names import RandomSampleNames
from random_tumor_types import RandomTumorTypes
import plot_data
import pandas as pd
import numpy as np
from scipy import stats
from multiprocessing import Pool
import itertools as it
import time
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
                                          rrclf.tsg_precision,
                                          rrclf.driver_precision],
                            'recall': [rrclf.onco_recall,
                                       rrclf.tsg_recall,
                                       rrclf.driver_recall],
                            'ROC AUC': [rrclf.onco_mean_roc_auc,
                                        rrclf.tsg_mean_roc_auc,
                                        rrclf.driver_mean_roc_auc],
                            'PR AUC': [rrclf.onco_mean_pr_auc,
                                       rrclf.tsg_mean_pr_auc,
                                       rrclf.driver_mean_pr_auc],
                            'count': [rrclf.onco_gene_count,
                                      rrclf.tsg_gene_count,
                                      rrclf.cancer_gene_count]},
                           index=['oncogene', 'tsg', 'driver'])
    return results


def py_random_forest(df, opts):
    rclf = RandomForest(df,
                        total_iter=1,
                        ntrees=opts['ntrees'],
                        min_ct=opts['min_count'])
    rclf.kfold_validation()  # run random forest
    results = pd.DataFrame({'precision': [rclf.onco_precision,
                                          rclf.tsg_precision,
                                          rclf.driver_precision],
                            'recall': [rclf.onco_recall,
                                       rclf.tsg_recall,
                                       rclf.driver_recall],
                            'ROC AUC': [rclf.onco_mean_roc_auc,
                                        rclf.tsg_mean_roc_auc,
                                        rclf.driver_mean_roc_auc],
                            'PR AUC': [rclf.onco_mean_pr_auc,
                                       rclf.tsg_mean_pr_auc,
                                       rclf.driver_mean_pr_auc],
                            'count': [rclf.onco_gene_count,
                                      rclf.tsg_gene_count,
                                      rclf.cancer_gene_count]},
                           index=['oncogene', 'tsg', 'driver'])
    return results


def naive_bayes(df, opts):
    nbclf = MultinomialNaiveBayes(df, min_ct=opts['min_count'], total_iter=1)
    nbclf.kfold_validation()  # run random forest
    results = pd.DataFrame({'precision': [nbclf.onco_precision,
                                          nbclf.tsg_precision,
                                          nbclf.driver_precision],
                            'recall': [nbclf.onco_recall,
                                       nbclf.tsg_recall,
                                       nbclf.driver_recall],
                            'ROC AUC': [nbclf.onco_mean_roc_auc,
                                        nbclf.tsg_mean_roc_auc,
                                        nbclf.driver_mean_roc_auc],
                            'PR AUC': [nbclf.onco_mean_pr_auc,
                                       nbclf.tsg_mean_pr_auc,
                                       nbclf.driver_mean_pr_auc],
                            'count': [nbclf.onco_gene_count,
                                      nbclf.tsg_gene_count,
                                      nbclf.cancer_gene_count]},
                           index=['oncogene', 'tsg', 'driver'])
    return results


def vogelstein_classification(df, opts):
    x, y_true = features.randomize(df)  # don't need randomization but gets x, y
    y_true_onco = (y_true==1).astype(int)
    y_true_tsg = (y_true==2).astype(int)
    y_true_driver = y_true_onco + y_true_tsg

    # input for 20/20 rule classifier
    ct_list = x[['recurrent count', 'deleterious count', 'total']].values

    # predict
    vclf = VogelsteinClassifier(min_count=opts['min_count'])
    y_pred = pd.Series(vclf.predict_list(ct_list,
                                         scale_type=opts['scale_type']))
    y_pred_onco = (y_pred=='oncogene').astype(int)
    y_pred_tsg = (y_pred=='tsg').astype(int)
    y_pred_driver = y_pred_onco + y_pred_tsg
    onco_prec, onco_recall, onco_fscore, onco_support = metrics.precision_recall_fscore_support(y_true_onco, y_pred_onco)
    tsg_prec, tsg_recall, tsg_fscore, tsg_support = metrics.precision_recall_fscore_support(y_true_tsg, y_pred_tsg)
    driver_prec, driver_recall, driver_fscore, driver_support = metrics.precision_recall_fscore_support(y_true_driver, y_pred_driver)

    # aggregate results
    results = pd.DataFrame({'precision': [onco_prec[1],
                                          tsg_prec[1],
                                          driver_prec[1]],
                            'recall': [onco_recall[1],
                                       tsg_recall[1],
                                       driver_recall[1]],
                            'count': [np.sum(y_pred_onco),
                                      np.sum(y_pred_tsg),
                                      np.sum(y_pred_driver)]},
                           index=['oncogene', 'tsg', 'driver'])

    return results


def retrieve_gene_features(opts):
    """Wrapper arround the retrieve_gene_features function in the
    features module."""
    # get additional features
    db_cfg = _utils.get_db_config('2020plus')
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
    return all_features


def calculate_sem(wp):
    """Calculates the standard error of the mean for a pd.Panel object.

    **Note:** The pd.Panel.apply method seems to have a bug preventing
    me from using it. So instead I am using the numpy apply function
    for calculating sem.

    Parameters
    ----------
    wp : pd.Panel
        panel that stratifies samples

    Returns
    -------
    tmp_sem : pd.DataFrame
        standard error of the mean calculated along the sample axis
    """
    tmp_sem_matrix = np.apply_along_axis(stats.sem, 0, wp.values)  # hack because pandas apply method has a bug
    tmp_sem = pd.DataFrame(tmp_sem_matrix,
                           columns=wp.minor_axis,
                           index=['oncogene', 'tsg', 'driver'])
    return tmp_sem


def calculate_stats(result_dict,
                    metrics=['precision', 'recall', 'ROC AUC', 'PR AUC', 'count']):
    """Computes mean and sem of classification performance metrics.

    Parameters
    ----------
    result_dict : dict
        dictionary with the i'th sample as the key and data frames
        with "oncogene"/"tsg" (row) classification performance metrics
        (columns) as values

    Returns
    -------
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


def save_simulation_result(r_panel, r_path,
                           py_panel, py_path,
                           nb_panel, nb_path,
                           v_panel, v_path,
                           v_lin_panel, v_lin_path):
    # make 'oncogenes'/'tsg' as the 'items' axis
    r_pan = r_panel.swapaxes('items', 'major')
    py_pan = py_panel.swapaxes('items', 'major')
    nb_pan = nb_panel.swapaxes('items', 'major')
    v_pan = v_panel.swapaxes('items', 'major')
    v_lin_pan = v_lin_panel.swapaxes('items', 'major')

    # collapse pd.panel into a data frame
    r_df = pd.merge(r_pan['oncogene'], r_pan['tsg'],
                    left_index=True, right_index=True,
                    suffixes=(' oncogene', ' tsg'))
    r_df = pd.merge(r_df, r_pan['driver'],
                    left_index=True, right_index=True,
                    suffixes=('', ' driver'))
    py_df = pd.merge(py_pan['oncogene'], py_pan['tsg'],
                     left_index=True, right_index=True,
                     suffixes=(' oncogene', ' tsg'))
    py_df = pd.merge(py_df, py_pan['driver'],
                    left_index=True, right_index=True,
                    suffixes=('', ' driver'))
    nb_df = pd.merge(nb_pan['oncogene'], nb_pan['tsg'],
                     left_index=True, right_index=True,
                     suffixes=(' oncogene', ' tsg'))
    nb_df = pd.merge(nb_df, nb_pan['driver'],
                    left_index=True, right_index=True,
                    suffixes=('', ' driver'))
    v_df = pd.merge(v_pan['oncogene'], v_pan['tsg'],
                    left_index=True, right_index=True,
                    suffixes=(' oncogene', ' tsg'))
    v_df = pd.merge(v_df, v_pan['driver'],
                    left_index=True, right_index=True,
                    suffixes=('', ' driver'))
    v_lin_df = pd.merge(v_lin_pan['oncogene'], v_lin_pan['tsg'],
                    left_index=True, right_index=True,
                    suffixes=(' oncogene', ' tsg'))
    v_lin_df = pd.merge(v_lin_df, v_lin_pan['driver'],
                    left_index=True, right_index=True,
                    suffixes=('', ' driver'))

    # save data frame to specified paths
    r_df.to_csv(r_path, sep='\t')
    py_df.to_csv(py_path, sep='\t')
    nb_df.to_csv(nb_path, sep='\t')
    v_df.to_csv(v_path, sep='\t')
    v_lin_df.to_csv(v_lin_path, sep='\t')


def predict(df, gene_df, opts):
    """Function called by multiprocessing to run predictions.

    **Parameters**
    info : tuple, length 3
        Elements:
          1. df, data frame of actual mutation data
          2. gene_df, data frame of gene features
          3. opts, dictionary of cli options from argparse
    """
    # df, gene_df, opts = info  # unpack tuple
    # df, opts = info  # unpack tuple

    #if opts['bootstrap']:
        # create feature matrix for bootstrap sample
        #df = features.process_features(df, 0)

    # merge gene features
    df = merge_feature_df(df, gene_df)  # all feature info
    df = df.set_index('gene')  # get rid of gene column

    # run classifiers on bootstrap sampled counts
    r_sim_results = r_random_forest(df, opts)
    py_sim_results = py_random_forest(df, opts)
    nb_sim_results = naive_bayes(df, opts)

    # do both non-scaling and linear scaling
    opts['scale_type'] = None
    v_sim_results = vogelstein_classification(df, opts)
    opts['scale_type'] = 'linear'
    v_linear_sim_results = vogelstein_classification(df, opts)

    return r_sim_results, py_sim_results, nb_sim_results, v_sim_results, v_linear_sim_results


def yield_sim_df(dfg, num_processes):
    """Yields (a generator) data frame objects according to
    the number of processes.

    **Parameters**

    dfg : one of my random sampling classes
        either a bootstrap object or RandomSampleNames object
    num_processes : int
        number of processes indicates how many data frames
        to yield.

    **returns**

    dfg_list : list
        a list of data frames of length num_processes
    """
    # yield lists of data frames until the last modulo
    # operator returns 0
    dfg_list = []
    for i, df in enumerate(dfg.dataframe_generator()):
        dfg_list.append(df)
        i_plus_1 = i + 1  # operator precedence requires this statement
        if not i_plus_1 % num_processes:
            yield dfg_list
            dfg_list = []

    # if number of data frames not perfectly divisible
    # by number of processes
    if dfg_list:
        yield dfg_list


def multiprocess_predict(dfg, gene_df, opts):
    num_processes = opts['processes']
    process_results = None

    for i in range(0, dfg.num_iter, num_processes):
        pool = Pool(processes=num_processes)
        del process_results  # possibly help free up more memory
        time.sleep(5)  # wait 5 seconds, might help make sure memory is free
        tmp_num_pred = dfg.num_iter - i if  i + num_processes > dfg.num_iter else num_processes
        # df_generator = dfg.dataframe_generator()
        info_repeat = it.repeat((dfg, gene_df, opts), tmp_num_pred)
        #pool = Pool(processes=tmp_num_pred)
        process_results = pool.imap(singleprocess_predict, info_repeat)
        pool.close()
        pool.join()
        yield process_results


def singleprocess_predict(info):
    dfg, gene_df, opts = info  # unpack tuple
    df = next(dfg.dataframe_generator())
    single_result = predict(df, gene_df, opts)
    return single_result


def main(cli_opts):
    out_opts = _utils.get_output_config('feature_matrix')
    sim_opts = _utils.get_output_config('simulation')
    my_db_path = _utils.get_db_config('2020plus')['db']
    conn = sqlite3.connect(my_db_path)

    gene_df = retrieve_gene_features(cli_opts)  # features like gene length, etc

    df = pd.read_csv(_utils.result_dir + out_opts['gene_feature_matrix'],
                     sep='\t', index_col=0)

    # iterate through each sampling rate
    r_result, py_result, nb_result, v_result, v_linear_result = {}, {}, {}, {}, {}
    for sample_rate in np.linspace(cli_opts['lower_sample_rate'],
                                   cli_opts['upper_sample_rate'],
                                   cli_opts['num_sample_rate']):
        if cli_opts['bootstrap']:
            # bootstrap object for sub-sampling
            dfg = Bootstrap(df.copy(),
                            subsample=sample_rate,
                            num_iter=cli_opts['samples'])
        elif cli_opts['random_samples']:
            # sample with out replacement of sample names
            dfg = RandomSampleNames(sub_sample=sample_rate,
                                    num_iter=cli_opts['samples'],
                                    db_conn=conn)
        else:
            # sample with out replacement of tumor types
            dfg = RandomTumorTypes(sub_sample=sample_rate,
                                   num_iter=cli_opts['samples'],
                                   db_conn=conn)


        # iterate through each sampled data set
        r_sim_results, py_sim_results, nb_sim_results, v_sim_results, v_linear_sim_results = {}, {}, {}, {}, {}
        i = 0  # counter for index of data frames
        for df_list in multiprocess_predict(dfg, gene_df, cli_opts):
            # save all the results
            for j, mydf in enumerate(df_list):
                r_sim_results[i],  py_sim_results[i], nb_sim_results[i], v_sim_results[i], v_linear_sim_results[i] = mydf
                i += 1

        # record result for a specific sample rate
        tmp_results = calculate_stats(r_sim_results)
        r_result[sample_rate] = tmp_results
        tmp_results = calculate_stats(py_sim_results)
        py_result[sample_rate] = tmp_results
        tmp_results = calculate_stats(nb_sim_results)
        nb_result[sample_rate] = tmp_results
        tmp_results = calculate_stats(v_sim_results, ['precision', 'recall', 'count'])
        v_result[sample_rate] = tmp_results
        tmp_results = calculate_stats(v_linear_sim_results, ['precision', 'recall', 'count'])
        v_linear_result[sample_rate] = tmp_results

    # make pandas panel objects out of summarized
    # results from simulations
    r_panel_result = pd.Panel(r_result)
    py_panel_result = pd.Panel(py_result)
    nb_panel_result = pd.Panel(nb_result)
    v_panel_result = pd.Panel(v_result)
    v_linear_panel_result = pd.Panel(v_linear_result)

    # save results of simulations in a text file
    r_path = _utils.sim_result_dir + sim_opts['r_result']
    py_path = _utils.sim_result_dir + sim_opts['py_result']
    nb_path = _utils.sim_result_dir + sim_opts['nb_result']
    v_path = _utils.sim_result_dir + sim_opts['v_result']
    v_linear_path = _utils.sim_result_dir + sim_opts['v_linear_result']
    save_simulation_result(r_panel_result, r_path,
                           py_panel_result, py_path,
                           nb_panel_result, nb_path,
                           v_panel_result, v_path,
                           v_linear_panel_result, v_linear_path)

    # aggregate results for plotting
    results = {'20/20+ classifier': r_panel_result}
               #'random forest': py_panel_result,
               #'naive bayes': nb_panel_result}

    # plot oncogene results of simulations
    tmp_save_path = _utils.sim_plot_dir + sim_opts['onco_pr_plot']
    plot_data.pr_auc_errorbar(results,
                              gene_type='oncogene',
                              save_path=tmp_save_path,
                              title='Oncogene PR AUC vs. DB size',
                              xlabel='Sample rate',
                              ylabel='PR AUC')
    tmp_save_path = _utils.sim_plot_dir + sim_opts['onco_roc_plot']
    plot_data.roc_auc_errorbar(results,
                               gene_type='oncogene',
                               save_path=tmp_save_path,
                               title='Oncogene ROC AUC vs. DB size',
                               xlabel='Sample rate',
                               ylabel='ROC AUC')

    # plot TSG results of simulations
    tmp_save_path = _utils.sim_plot_dir + sim_opts['tsg_pr_plot']
    plot_data.pr_auc_errorbar(results,
                              gene_type='tsg',
                              save_path=tmp_save_path,
                              title='TSG PR AUC vs. DB size',
                              xlabel='Sample rate',
                              ylabel='PR AUC')
    tmp_save_path = _utils.sim_plot_dir + sim_opts['tsg_roc_plot']
    plot_data.roc_auc_errorbar(results,
                               gene_type='tsg',
                               save_path=tmp_save_path,
                               title='TSG ROC AUC vs. DB size',
                               xlabel='Sample rate',
                               ylabel='ROC AUC')

    # plot driver auc metrics
    tmp_save_path = _utils.sim_plot_dir + sim_opts['driver_pr_plot']
    plot_data.pr_auc_errorbar(results,
                              gene_type='driver',
                              save_path=tmp_save_path,
                              title='Driver PR AUC vs. DB size',
                              xlabel='Sample rate',
                              ylabel='PR AUC')
    tmp_save_path = _utils.sim_plot_dir + sim_opts['driver_roc_plot']
    plot_data.roc_auc_errorbar(results,
                               gene_type='driver',
                               save_path=tmp_save_path,
                               title='Driver ROC AUC vs. DB size',
                               xlabel='Sample rate',
                               ylabel='ROC AUC')

    # since the vogelstein classifier doesn't predict probabilities
    # I can't generate a PR or ROC curve. However, I can evaluate
    # metrics like precision and recall
    # results['20/20 rule'] = v_panel_result
    # results['20/20 rule (linear scaling)'] = v_linear_panel_result

    # plot oncogene precision/recall
    tmp_save_path = _utils.sim_plot_dir + sim_opts['onco_precision_plot']
    plot_data.precision_errorbar(results,
                                 gene_type='oncogene',
                                 save_path=tmp_save_path,
                                 title='Oncogene Precision vs. DB size',
                                 xlabel='Sample rate',
                                 ylabel='Precision')
    tmp_save_path = _utils.sim_plot_dir + sim_opts['onco_recall_plot']
    plot_data.recall_errorbar(results,
                              gene_type='oncogene',
                              save_path=tmp_save_path,
                              title='Oncogene Recall vs. DB size',
                              xlabel='Sample rate',
                              ylabel='Recall')

    # plot TSG precision/recall
    tmp_save_path = _utils.sim_plot_dir + sim_opts['tsg_precision_plot']
    plot_data.precision_errorbar(results,
                                 gene_type='tsg',
                                 save_path=tmp_save_path,
                                 title='TSG Precision vs. DB size',
                                 xlabel='Sample rate',
                                 ylabel='Precision')
    tmp_save_path = _utils.sim_plot_dir + sim_opts['tsg_recall_plot']
    plot_data.recall_errorbar(results,
                              gene_type='tsg',
                              save_path=tmp_save_path,
                              title='TSG Recall vs. DB size',
                              xlabel='Sample rate',
                              ylabel='Recall')

    # plot driver precision/recall
    tmp_save_path = _utils.sim_plot_dir + sim_opts['driver_precision_plot']
    plot_data.precision_errorbar(results,
                                 gene_type='driver',
                                 save_path=tmp_save_path,
                                 title='Driver Precision vs. DB size',
                                 xlabel='Sample rate',
                                 ylabel='Precision')
    tmp_save_path = _utils.sim_plot_dir + sim_opts['driver_recall_plot']
    plot_data.recall_errorbar(results,
                              gene_type='driver',
                              save_path=tmp_save_path,
                              title='Driver Recall vs. DB size',
                              xlabel='Sample rate',
                              ylabel='Recall')

    # delete naive bayes since it predicts a lot of genes
    # del results['naive bayes']
    # del results['random forest']

    # plot number of predicted genes
    tmp_save_path = _utils.sim_plot_dir + sim_opts['onco_count_plot']
    plot_data.count_errorbar(results,
                             gene_type='oncogene',
                             save_path=tmp_save_path,
                             title='Number of predicted oncogenes vs. DB size',
                             xlabel='Sample rate',
                             ylabel='Number of predicted oncogenes')
    tmp_save_path = _utils.sim_plot_dir + sim_opts['tsg_count_plot']
    plot_data.count_errorbar(results,
                             gene_type='tsg',
                             save_path=tmp_save_path,
                             title='Number of predicted TSG vs. DB size',
                             xlabel='Sample rate',
                             ylabel='Number of predicted TSG')
    tmp_save_path = _utils.sim_plot_dir + sim_opts['driver_count_plot']
    plot_data.count_errorbar(results,
                             gene_type='driver',
                             save_path=tmp_save_path,
                             title='Number of predicted driver vs. DB size',
                             xlabel='Sample rate',
                             ylabel='Number of predicted driver')
