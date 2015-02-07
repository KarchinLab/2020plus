from src.classify.python.r_random_forest_clf import RRandomForest
from src.classify.python.random_forest_clf import RandomForest
from src.classify.python.multinomial_nb_clf import MultinomialNaiveBayes
import src.simulation.python.simulation as sim
import src.utils.python.util as _utils
import simulate_performance as sp
from random_split import RandomSplit
import pandas as pd
import scipy.stats as stats
import sqlite3
from multiprocessing import Pool
import itertools as it
import time


def r_random_forest_compare(df1, df2, opts):
    # predict and compile results
    rrclf1 = RRandomForest(df1,
                           other_sample_ratio=opts['other_ratio'],
                           driver_sample=opts['driver_rate'],
                           ntrees=opts['ntrees'],
                           min_ct=opts['min_count'])
    onco_prob1, tsg_prob1, other_prob1 = rrclf1.kfold_prediction()
    driver_prob1 = onco_prob1 + tsg_prob1
    df1 = pd.DataFrame({'onco': onco_prob1,
                        'tsg': tsg_prob1,
                        'other': other_prob1,
                        'driver': driver_prob1},
                       index=onco_prob1.index)
    rrclf2 = RRandomForest(df2,
                           other_sample_ratio=opts['other_ratio'],
                           driver_sample=opts['driver_rate'],
                           ntrees=opts['ntrees'],
                           min_ct=opts['min_count'])
    onco_prob2, tsg_prob2, other_prob2 = rrclf2.kfold_prediction()
    driver_prob2 = onco_prob2 + tsg_prob2
    df2 = pd.DataFrame({'onco': onco_prob2,
                        'tsg': tsg_prob2,
                        'other': other_prob2,
                        'driver': driver_prob2},
                       index=onco_prob2.index)
    # get predicted class
    onco_pred1 = (df1['onco'] > df1['tsg']) & (df1['onco'] > df1['other'])
    onco_pred2 = (df2['onco'] > df2['tsg']) & (df2['onco'] > df2['other'])
    tsg_pred1 = (df1['tsg'] > df1['onco']) & (df1['tsg'] > df1['other'])
    tsg_pred2 = (df2['tsg'] > df2['onco']) & (df2['tsg'] > df2['other'])

    # rank genes
    driver_top_rank1, driver_top_rank2 = sim.rank_genes(driver_prob1, driver_prob2, thresh=.5)
    onco_top_rank1, onco_top_rank2 = sim.rank_genes(onco_prob1, onco_prob2,
                                                    mask1=onco_pred1, mask2=onco_pred2)
    tsg_top_rank1, tsg_top_rank2 = sim.rank_genes(tsg_prob1, tsg_prob2,
                                                  mask1=tsg_pred1, mask2=tsg_pred2)

    # calculate jaccard index
    driver_jaccard = sim.overlap(driver_prob1, driver_prob2, thresh=.5)
    onco_jaccard = sim.overlap(onco_prob1, onco_prob2,
                                     mask1=onco_pred1, mask2=onco_pred2)
    tsg_jaccard = sim.overlap(tsg_prob1, tsg_prob2,
                                    mask1=tsg_pred1, mask2=tsg_pred2)

    # calc spearman rank correlation
    sp_rho_driver, sp_pval_driver = stats.pearsonr(driver_top_rank1, driver_top_rank2)
    sp_rho_onco, sp_pval_onco = stats.pearsonr(onco_top_rank1, onco_top_rank2)
    sp_rho_tsg, sp_pval_tsg = stats.pearsonr(tsg_top_rank1, tsg_top_rank2)

    # calc kendall tau correlation
    #kt_rho_driver, kt_pval_driver = stats.kendalltau(driver_top_rank1, driver_top_rank2)
    #kt_rho_onco, kt_pval_onco = stats.kendalltau(onco_top_rank1, onco_top_rank2)
    #kt_rho_tsg, kt_pval_tsg = stats.kendalltau(tsg_top_rank1, tsg_top_rank2)

    # calculate jaccard index at specified depths
    driver_prob1.sort(ascending=False)
    driver_prob2.sort(ascending=False)
    tmp_tuple = sim.weighted_overlap(driver_prob1, driver_prob2,
                                           max_depth=opts['depth'],
                                           step_size=opts['step_size'],
                                           weight_factor=opts['weight'])
    _, mean_driver_ji, weighted_mean_driver_ji = tmp_tuple
    onco_prob1.sort(ascending=False)
    onco_prob2.sort(ascending=False)
    tmp_tuple = sim.weighted_overlap(onco_prob1, onco_prob2,
                                           max_depth=opts['depth'],
                                           step_size=opts['step_size'],
                                           weight_factor=opts['weight'])
    _, mean_onco_ji, weighted_mean_onco_ji = tmp_tuple
    tsg_prob1.sort(ascending=False)
    tsg_prob2.sort(ascending=False)
    tmp_tuple = sim.weighted_overlap(tsg_prob1, tsg_prob2,
                                           max_depth=opts['depth'],
                                           step_size=opts['step_size'],
                                           weight_factor=opts['weight'])
    _, mean_tsg_ji, weighted_mean_tsg_ji = tmp_tuple

    results = pd.DataFrame({'jaccard index': [onco_jaccard,
                                              tsg_jaccard,
                                              driver_jaccard],
                            'spearman correlation': [sp_rho_onco,
                                                     sp_rho_tsg,
                                                     sp_rho_driver],
                            #'kendall tau correlation': [kt_rho_onco,
                                                        #kt_rho_tsg,
                                                        #kt_rho_driver],
                            'spearman p-value': [sp_pval_onco,
                                                 sp_pval_tsg,
                                                 sp_pval_driver],
                            #'kendall tau p-value': [kt_pval_onco,
                                                    #kt_pval_tsg,
                                                    #kt_pval_driver],
                            'mean jaccard index': [mean_onco_ji,
                                                   mean_tsg_ji,
                                                   mean_driver_ji],
                            'weighted mean jaccard index': [weighted_mean_onco_ji,
                                                            weighted_mean_tsg_ji,
                                                            weighted_mean_driver_ji]},
                            index=['oncogene', 'tsg', 'driver'])
    return results


def py_random_forest_compare(df1, df2, opts):
    rclf1 = RandomForest(df1,
                         ntrees=opts['ntrees'],
                         min_ct=opts['min_count'])
    onco_prob1, tsg_prob1, other_prob1 = rclf1.kfold_prediction()
    driver_prob1 = onco_prob1 + tsg_prob1
    df1 = pd.DataFrame({'onco': onco_prob1,
                        'tsg': tsg_prob1,
                        'other': other_prob1,
                        'driver': driver_prob1},
                       index=onco_prob1.index)
    rclf2 = RandomForest(df2,
                         ntrees=opts['ntrees'],
                         min_ct=opts['min_count'])
    onco_prob2, tsg_prob2, other_prob2 = rclf2.kfold_prediction()
    driver_prob2 = onco_prob2 + tsg_prob2
    df2 = pd.DataFrame({'onco': onco_prob2,
                        'tsg': tsg_prob2,
                        'other': other_prob2,
                        'driver': driver_prob2},
                       index=onco_prob2.index)

    # get predicted class
    onco_pred1 = (df1['onco'] > df1['tsg']) & (df1['onco'] > df1['other'])
    onco_pred2 = (df2['onco'] > df2['tsg']) & (df2['onco'] > df2['other'])
    tsg_pred1 = (df1['tsg'] > df1['onco']) & (df1['tsg'] > df1['other'])
    tsg_pred2 = (df2['tsg'] > df2['onco']) & (df2['tsg'] > df2['other'])

    # rank genes
    driver_top_rank1, driver_top_rank2 = sim.rank_genes(driver_prob1, driver_prob2, thresh=.5)
    onco_top_rank1, onco_top_rank2 = sim.rank_genes(onco_prob1, onco_prob2,
                                                    mask1=onco_pred1, mask2=onco_pred2)
    tsg_top_rank1, tsg_top_rank2 = sim.rank_genes(tsg_prob1, tsg_prob2,
                                                  mask1=tsg_pred1, mask2=tsg_pred2)

    # calculate jaccard index
    driver_jaccard = sim.overlap(driver_prob1, driver_prob2, thresh=.5)
    onco_jaccard = sim.overlap(onco_prob1, onco_prob2,
                                     mask1=onco_pred1, mask2=onco_pred2)
    tsg_jaccard = sim.overlap(tsg_prob1, tsg_prob2,
                                    mask1=tsg_pred1, mask2=tsg_pred2)

    # calc spearman rank correlation
    sp_rho_driver, sp_pval_driver = stats.pearsonr(driver_top_rank1, driver_top_rank2)
    sp_rho_onco, sp_pval_onco = stats.pearsonr(onco_top_rank1, onco_top_rank2)
    sp_rho_tsg, sp_pval_tsg = stats.pearsonr(tsg_top_rank1, tsg_top_rank2)

    # calc kendall tau correlation
    #kt_rho_driver, kt_pval_driver = stats.kendalltau(driver_top_rank1, driver_top_rank2)
    #kt_rho_onco, kt_pval_onco = stats.kendalltau(onco_top_rank1, onco_top_rank2)
    #kt_rho_tsg, kt_pval_tsg = stats.kendalltau(tsg_top_rank1, tsg_top_rank2)

    # calculate jaccard index at specified depths
    driver_prob1.sort(ascending=False)
    driver_prob2.sort(ascending=False)
    tmp_tuple = sim.weighted_overlap(driver_prob1, driver_prob2,
                                           max_depth=opts['depth'],
                                           step_size=opts['step_size'],
                                           weight_factor=opts['weight'])
    _, mean_driver_ji, weighted_mean_driver_ji = tmp_tuple
    onco_prob1.sort(ascending=False)
    onco_prob2.sort(ascending=False)
    tmp_tuple = sim.weighted_overlap(onco_prob1, onco_prob2,
                                           max_depth=opts['depth'],
                                           step_size=opts['step_size'],
                                           weight_factor=opts['weight'])
    _, mean_onco_ji, weighted_mean_onco_ji = tmp_tuple
    tsg_prob1.sort(ascending=False)
    tsg_prob2.sort(ascending=False)
    tmp_tuple = sim.weighted_overlap(tsg_prob1, tsg_prob2,
                                           max_depth=opts['depth'],
                                           step_size=opts['step_size'],
                                           weight_factor=opts['weight'])
    _, mean_tsg_ji, weighted_mean_tsg_ji = tmp_tuple

    results = pd.DataFrame({'jaccard index': [onco_jaccard,
                                              tsg_jaccard,
                                              driver_jaccard],
                            'spearman correlation': [sp_rho_onco,
                                                     sp_rho_tsg,
                                                     sp_rho_driver],
                            #'kendall tau correlation': [kt_rho_onco,
                                                        #kt_rho_tsg,
                                                        #kt_rho_driver],
                            'spearman p-value': [sp_pval_onco,
                                                 sp_pval_tsg,
                                                 sp_pval_driver],
                            #'kendall tau p-value': [kt_pval_onco,
                                                    #kt_pval_tsg,
                                                    #kt_pval_driver],
                            'mean jaccard index': [mean_onco_ji,
                                                   mean_tsg_ji,
                                                   mean_driver_ji],
                            'weighted mean jaccard index': [weighted_mean_onco_ji,
                                                            weighted_mean_tsg_ji,
                                                            weighted_mean_driver_ji]},
                            index=['oncogene', 'tsg', 'driver'])
    return results


def naive_bayes_compare(df1, df2, opts):
    nbclf1 = MultinomialNaiveBayes(df1, min_ct=opts['min_count'])
    onco_prob1, tsg_prob1, other_prob1 = nbclf1.kfold_prediction()
    driver_prob1 = onco_prob1 + tsg_prob1
    df1 = pd.DataFrame({'onco': onco_prob1,
                        'tsg': tsg_prob1,
                        'other': other_prob1,
                        'driver': driver_prob1},
                       index=onco_prob1.index)
    nbclf2 = MultinomialNaiveBayes(df2, min_ct=opts['min_count'])
    onco_prob2, tsg_prob2, other_prob2 = nbclf2.kfold_prediction()
    driver_prob2 = onco_prob2 + tsg_prob2
    df2 = pd.DataFrame({'onco': onco_prob2,
                        'tsg': tsg_prob2,
                        'other': other_prob2,
                        'driver': driver_prob2},
                       index=onco_prob2.index)

    # get predicted class
    onco_pred1 = (df1['onco'] > df1['tsg']) & (df1['onco'] > df1['other'])
    onco_pred2 = (df2['onco'] > df2['tsg']) & (df2['onco'] > df2['other'])
    tsg_pred1 = (df1['tsg'] > df1['onco']) & (df1['tsg'] > df1['other'])
    tsg_pred2 = (df2['tsg'] > df2['onco']) & (df2['tsg'] > df2['other'])

    # rank genes
    driver_top_rank1, driver_top_rank2 = sim.rank_genes(driver_prob1, driver_prob2, thresh=.5)
    onco_top_rank1, onco_top_rank2 = sim.rank_genes(onco_prob1, onco_prob2,
                                                    mask1=onco_pred1, mask2=onco_pred2)
    tsg_top_rank1, tsg_top_rank2 = sim.rank_genes(tsg_prob1, tsg_prob2,
                                                  mask1=tsg_pred1, mask2=tsg_pred2)

    # calculate jaccard index
    driver_jaccard = sim.overlap(driver_prob1, driver_prob2, thresh=.5)
    onco_jaccard = sim.overlap(onco_prob1, onco_prob2,
                                     mask1=onco_pred1, mask2=onco_pred2)
    tsg_jaccard = sim.overlap(tsg_prob1, tsg_prob2,
                                    mask1=tsg_pred1, mask2=tsg_pred2)

    # calc spearman rank correlation
    sp_rho_driver, sp_pval_driver = stats.pearsonr(driver_top_rank1, driver_top_rank2)
    sp_rho_onco, sp_pval_onco = stats.pearsonr(onco_top_rank1, onco_top_rank2)
    sp_rho_tsg, sp_pval_tsg = stats.pearsonr(tsg_top_rank1, tsg_top_rank2)

    # calc kendall tau correlation
    #kt_rho_driver, kt_pval_driver = stats.kendalltau(driver_top_rank1, driver_top_rank2)
    #kt_rho_onco, kt_pval_onco = stats.kendalltau(onco_top_rank1, onco_top_rank2)
    #kt_rho_tsg, kt_pval_tsg = stats.kendalltau(tsg_top_rank1, tsg_top_rank2)

    # calculate jaccard index at specified depths
    driver_prob1.sort(ascending=False)
    driver_prob2.sort(ascending=False)
    tmp_tuple = sim.weighted_overlap(driver_prob1, driver_prob2,
                                           max_depth=opts['depth'],
                                           step_size=opts['step_size'],
                                           weight_factor=opts['weight'])
    _, mean_driver_ji, weighted_mean_driver_ji = tmp_tuple
    onco_prob1.sort(ascending=False)
    onco_prob2.sort(ascending=False)
    tmp_tuple = sim.weighted_overlap(onco_prob1, onco_prob2,
                                           max_depth=opts['depth'],
                                           step_size=opts['step_size'],
                                           weight_factor=opts['weight'])
    _, mean_onco_ji, weighted_mean_onco_ji = tmp_tuple
    tsg_prob1.sort(ascending=False)
    tsg_prob2.sort(ascending=False)
    tmp_tuple = sim.weighted_overlap(tsg_prob1, tsg_prob2,
                                           max_depth=opts['depth'],
                                           step_size=opts['step_size'],
                                           weight_factor=opts['weight'])
    _, mean_tsg_ji, weighted_mean_tsg_ji = tmp_tuple

    results = pd.DataFrame({'jaccard index': [onco_jaccard,
                                              tsg_jaccard,
                                              driver_jaccard],
                            'spearman correlation': [sp_rho_onco,
                                                     sp_rho_tsg,
                                                     sp_rho_driver],
                            #'kendall tau correlation': [kt_rho_onco,
                                                        #kt_rho_tsg,
                                                        #kt_rho_driver],
                            'spearman p-value': [sp_pval_onco,
                                                 sp_pval_tsg,
                                                 sp_pval_driver],
                            #'kendall tau p-value': [kt_pval_onco,
                                                    #kt_pval_tsg,
                                                    #kt_pval_driver],
                            'mean jaccard index': [mean_onco_ji,
                                                   mean_tsg_ji,
                                                   mean_driver_ji],
                            'weighted mean jaccard index': [weighted_mean_onco_ji,
                                                            weighted_mean_tsg_ji,
                                                            weighted_mean_driver_ji]},
                            index=['oncogene', 'tsg', 'driver'])
    return results


def compare_prediction(df1, df2, gene_df, opts):
    """Function called by multiprocessing to run predictions and compare
    their predictions.

    Parameters
    ----------
    df1 : pd.DataFrame
        first half of the split of the gene features
    df2 : pd.DataFrame
        second half of the split of the gene features
    gene_df : pd.DataFrame
        not mutational features related to genes (e.g. gene length)
    opts : dict
        parameter options for classifiers

    Returns
    -------
    tuple
        results of predictions
    """
    # merge gene features
    df1 = pd.merge(df1, gene_df,
                   how='left', on='gene')
    df1 = df1.set_index('gene')
    df2 = pd.merge(df2, gene_df,
                   how='left', on='gene')
    df2 = df2.set_index('gene')

    # run classifiers on random split
    r_sim_results = r_random_forest_compare(df1, df2, opts)
    py_sim_results = py_random_forest_compare(df1, df2, opts)
    nb_sim_results = naive_bayes_compare(df1, df2, opts)

    return r_sim_results, py_sim_results, nb_sim_results


@_utils.log_error_decorator
def singleprocess_compare_prediction(info):
    dfg, gene_df, opts = info  # unpack tuple
    df1, df2 = next(dfg.dataframe_generator())
    single_result = compare_prediction(df1, df2, gene_df, opts)
    return single_result


def multiprocess_compare_prediction(dfg, gene_df, opts):
    # handle the number of processes to use, should it even use the
    # multiprocessing module?
    multiprocess_flag = opts['processes']>0
    if multiprocess_flag:
        num_processes = opts['processes']
    else:
        num_processes = 1
    opts['processes'] = 0  # do not use multi-processing within permutation test
    # num_processes = opts['processes']

    # initialize output
    process_results = None

    for i in range(0, dfg.num_iter, num_processes):
        if multiprocess_flag:
            pool = Pool(processes=num_processes)
            del process_results  # possibly help free up more memory
            time.sleep(5)  # wait 5 seconds, might help make sure memory is free
            tmp_num_pred = dfg.num_iter - i if  i + num_processes > dfg.num_iter else num_processes
            # df_generator = dfg.dataframe_generator()
            info_repeat = it.repeat((dfg, gene_df, opts), tmp_num_pred)
            #pool = Pool(processes=tmp_num_pred)
            process_results = pool.imap(singleprocess_compare_prediction, info_repeat)
            pool.close()
            pool.join()
            yield process_results
        else:
            info = (dfg, gene_df, opts)
            tmp_result = singleprocess_compare_prediction(info)
            yield (tmp_result,)


def main(cli_opts):
    out_opts = _utils.get_output_config('feature_matrix')
    sim_opts = _utils.get_output_config('simulate_consistency')
    my_db_path = _utils.get_db_config('2020plus')['db']
    conn = sqlite3.connect(my_db_path)

    gene_df = sp.retrieve_gene_features(cli_opts)  # features like gene length, etc

    r_result, py_result, nb_result = {}, {}, {}

    # object that generates features from randomly choosen sample names while
    # still respecting the stratification of tumor types
    if not cli_opts['with_replacement']:
        sample_rate = .5  # always do half splits
        dfg = RandomSplit(sub_sample=sample_rate,
                          num_iter=cli_opts['samples'],
                          db_conn=conn)
    else:
        sample_rate = cli_opts['with_replacement']
        dfg = RandomSplit(sub_sample=sample_rate,
                          num_iter=cli_opts['samples'],
                          db_conn=conn,
                          with_replacement=True)

    # iterate through each sampled data set
    r_sim_results, py_sim_results, nb_sim_results = {}, {}, {}
    i = 0  # counter for index of data frames
    for df_list in multiprocess_compare_prediction(dfg, gene_df, cli_opts):
        # save all the results
        for j, mydf in enumerate(df_list):
            r_sim_results[i],  py_sim_results[i], nb_sim_results[i] = mydf
            i += 1

    # record result for a specific sample rate
    tmp_results = sp.calculate_stats(r_sim_results)
    r_result[sample_rate] = tmp_results
    tmp_results = sp.calculate_stats(py_sim_results)
    py_result[sample_rate] = tmp_results
    tmp_results = sp.calculate_stats(nb_sim_results)
    nb_result[sample_rate] = tmp_results

    # make pandas panel objects out of summarized
    # results from simulations
    r_panel_result = pd.Panel(r_result)
    py_panel_result = pd.Panel(py_result)
    nb_panel_result = pd.Panel(nb_result)
    r_path = _utils.sim_result_dir + sim_opts['r_result']
    py_path = _utils.sim_result_dir + sim_opts['py_result']
    nb_path = _utils.sim_result_dir + sim_opts['nb_result']
    r_panel_result[sample_rate].to_csv(r_path, sep='\t')
    py_panel_result[sample_rate].to_csv(py_path, sep='\t')
    nb_panel_result[sample_rate].to_csv(nb_path, sep='\t')

    # results for plotting
    results = {'20/20+ classifier': r_panel_result[sample_rate],
               'Random Forest': py_panel_result[sample_rate],
               'Naive Bayes': nb_panel_result[sample_rate]}
