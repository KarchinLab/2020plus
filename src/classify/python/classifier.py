from __future__ import division
from vogelstein_classifier import VogelsteinClassifier
from random_forest_clf import RandomForest
from multinomial_nb_clf import MultinomialNaiveBayes
from dummy_clf import DummyClf
from r_random_forest_clf import RRandomForest
import src.utils.python.util as _utils
import plot_data
import pandas as pd
import numpy as np
import glob
import re
import bisect
import logging

logger = logging.getLogger(__name__)

def calc_class_info(df, onco_pct, tsg_pct, min_ct, tsg_min=None, kind='oncogene'):
    # calculate the number of genes classified as oncogene
    if not tsg_min:
        vclf = VogelsteinClassifier(onco_pct, tsg_pct, min_count=min_ct)
    else:
        vclf = VogelsteinClassifier(onco_pct, tsg_pct,
                                    min_count=min_ct,
                                    tsg_min=tsg_min)
    df['total'] = df.T.sum()
    input_list = [(row['recurrent missense'] + row['recurrent indel'],
                   row['frame shift'] + row['nonsense'] + row['lost stop'] + row['no protein'] + row['splicing mutation'],
                   row['total'])
                  for i, row in df.iterrows()]
    df['20/20 predicted class'] = vclf.predict_list(input_list)
    class_cts = df['20/20 predicted class'].value_counts()

    # calculate the pct of known oncogenes found
    df['true class'] = [_utils.classify_gene(gene)
                        for gene in df.index.tolist()]
    tmpdf = df.copy()  # prevent wierd behavior
    known_class = tmpdf[tmpdf['true class']==kind]
    num_class_found = len(known_class[known_class['20/20 predicted class']==kind])
    total_class = len(known_class)  # total number of oncogenes with counts
    pct_class_found = num_class_found / total_class

    return class_cts[kind], pct_class_found


def num_onco_by_recurrent_mutations(onco_pct, tsg_pct, min_ct):
    """Count number of oncogenes while varying the definition of recurrency"""
    # calculate counts for oncogenes/tsg with varying the required the number
    # of mutations to define a recurrent position
    file_match_pattern = _utils.result_dir + 'gene_feature_matrix.r*.txt'
    gene_design_matrix_paths = glob.glob(file_match_pattern)
    onco_ct_list, onco_pct_list = [], []  # list of cts/pct for oncogenes
    for file_path in gene_design_matrix_paths:
        tmp_df = pd.read_csv(file_path, sep='\t', index_col=0)
        tmp_ct, tmp_pct = calc_class_info(tmp_df,
                                         onco_pct=onco_pct,  # pct thresh for onco
                                         tsg_pct=tsg_pct,  # pct thresh for tsg
                                         min_ct=min_ct,  # min count for a gene
                                         kind='oncogene')
        onco_ct_list.append(tmp_ct)
        onco_pct_list.append(tmp_pct)

    # extract the '-r' parameter from the file name
    recur_param_pattern = '\d+'
    recur_param_list = [int(re.search(recur_param_pattern, mypath).group())
                        for mypath in gene_design_matrix_paths]

    # return dataframe with counts for each use of a recurrent mutation counts
    mycts = pd.Series(onco_ct_list, index=recur_param_list)
    mypct = pd.Series(onco_pct_list, index=recur_param_list)
    return mycts, mypct


def num_tsg_by_threshold(onco_pct, tsg_pct, min_ct, del_param_list):
    file_path = _utils.result_dir + 'gene_feature_matrix.txt'
    tmp_df = pd.read_csv(file_path, sep='\t', index_col=0)
    tsg_ct_list, tsg_pct_list = [], []
    for tmp_tsg_min in del_param_list:
        tmp_tsg_ct, tmp_tsg_pct = calc_class_info(tmp_df.copy(),
                                                  onco_pct=onco_pct,  # pct thresh for onco
                                                  tsg_pct=tsg_pct,  # pct thresh for tsg
                                                  min_ct=min_ct,  # min count for a gene
                                                  tsg_min=tmp_tsg_min,  # min del ct for tsg
                                                  kind='tsg')
        tsg_ct_list.append(tmp_tsg_ct)
        tsg_pct_list.append(tmp_tsg_pct)

    # return dataframe with counts for each use of a recurrent mutation counts
    mycts = pd.Series(tsg_ct_list, index=del_param_list)
    mypct = pd.Series(tsg_pct_list, index=del_param_list)
    return mycts, mypct


def num_pred_by_pct_threshold(min_ct):
    """Enumerates percent thresholds for 2020 rule."""
    # initialization of dataframe
    cts, pct = num_onco_by_recurrent_mutations(.2, .2, min_ct)
    onco_ct = pd.DataFrame(index=cts.index)
    onco_pct = pd.DataFrame(index=pct.index)

    # test different percentage thresholds
    thresh_range = np.arange(.15, .5, .05)
    tsg_ct_params = range(5, 17, 2)
    tsg_ct = pd.DataFrame(index=tsg_ct_params)
    tsg_pct = pd.DataFrame(index=tsg_ct_params)
    for threshold in thresh_range:
        tmp_ct, tmp_pct = num_onco_by_recurrent_mutations(threshold, threshold, min_ct)
        onco_ct[str(threshold)] = tmp_ct
        onco_pct[str(threshold)] = tmp_pct
        tmp_ct, tmp_pct = num_tsg_by_threshold(threshold, threshold,
                                               min_ct, tsg_ct_params)
        tsg_ct[threshold] = tmp_ct
        tsg_pct[threshold] = tmp_pct

    return onco_ct, onco_pct, tsg_ct, tsg_pct


def generate_2020_result(onco_pct, tsg_pct, min_ct):
    """Simply runs the 20/20 rule with given parameters."""
    # process count features from "data_analysis" results
    in_cfg = _utils.get_input_config('classifier')  # get directory
    df = pd.read_csv(_utils.save_dir + in_cfg['gene_feature'],
                     sep='\t', index_col=0)
    # df['total'] = df.T.sum()
    #df['total recurrent count'] = df['total'] * (df['recurrent missense'] + df['recurrent indel'])
    #df['total deleterious count'] = df['total'] * (df['frame shift'] + df['nonsense'] + df['lost stop'] + df['no protein'] + df['splicing mutation'])
    df['oncogene score'] = df['recurrent missense'] + df['recurrent indel']
    df['tsg score'] = df['frame shift'] + df['nonsense'] + df['lost stop'] + df['no protein'] + df['splicing mutation']
    #df['oncogene score'] = df['total recurrent count'].astype(float).div(df['total'])
    #df['tsg score'] = df['total deleterious count'].astype(float).div(df['total'])

    # predict using the "20/20" rule
    vclf = VogelsteinClassifier(onco_pct, tsg_pct, min_count=min_ct)
    input_list = [(row['recurrent count'], row['deleterious count'], row['total'])
                  for i, row in df.iterrows()]
    df['20/20 predicted class'] = vclf.predict_list(input_list)
    df['true class'] = [_utils.classify_gene(gene)
                        for gene in df.index.tolist()]
    return df


def compute_p_value(scores, empirical_p_values):
    """Get the p-value for each score by examining the list null distribution
    where scores are obtained by a certain probability.
    """
    num_scores = len(scores)
    pvals = pd.Series(np.zeros(num_scores))
    emp_p_val_scores = list(reversed(empirical_p_values.index.tolist()))
    num_emp_p_val_scores = len(emp_p_val_scores)
    score2pval = lambda x: empirical_p_values.iloc[num_emp_p_val_scores-max(bisect.bisect_right(emp_p_val_scores, x), 1)]
    pvals = scores.apply(score2pval)
    return pvals


def rand_forest_pred(clf, data, result_path, null_dist=None):
    """Makes gene predictions using a random forest classifier.

    Parameters
    ----------
    clf : GenericClassifier
        random forest sub-class of GenericClassifier
    data : pd.DataFrame
        data frame containing feature information
    result_path : str
        path to save text file result

    Returns
    -------
    tmp_df : pd.DataFrame
        random forest results (already saved to file)
    """
    # perform prediction
    onco_prob, tsg_prob, other_prob = clf.kfold_prediction()
    true_class = clf.y

    # save features/prediction results
    tmp_df = data.copy()
    tmp_df['oncogene score'] = onco_prob
    tmp_df['tsg score'] = tsg_prob
    tmp_df['other score'] = other_prob
    tmp_df['driver score'] = 1 - other_prob
    pred_class = tmp_df[['other score', 'oncogene score', 'tsg score']].values.argmax(axis=1)
    tmp_df['majority vote class'] = pred_class
    tmp_df['majority vote cancer gene'] = (tmp_df['driver score'] > .5).astype(int)
    tmp_df['training list class'] = true_class
    tmp_df = tmp_df.fillna(0)
    tmp_df = tmp_df.sort(['driver score',], ascending=False)

    if null_dist is not None:
        # add oncogene p-value
        onco_score = tmp_df['oncogene score'].copy()
        onco_score.sort(ascending=False)
        tmp_df['oncogene p-value'] = compute_p_value(onco_score,
                                                     null_dist['oncogene p-value'].dropna())
        tmp_df['oncogene q-value'] = _utils.bh_fdr(tmp_df['oncogene p-value'])

        # add tsg p-value
        tsg_score = tmp_df['tsg score'].copy()
        tsg_score.sort(ascending=False)
        tmp_df['tsg p-value'] = compute_p_value(tsg_score,
                                                null_dist['tsg p-value'].dropna())
        tmp_df['tsg q-value'] = _utils.bh_fdr(tmp_df['tsg p-value'])

        # add driver p-values
        tmp_df['driver p-value'] = compute_p_value(tmp_df['driver score'],
                                                   null_dist['driver p-value'].dropna())
        tmp_df['driver q-value'] = _utils.bh_fdr(tmp_df['driver p-value'])

    tmp_df.to_csv(result_path, sep='\t')

    return tmp_df


def trained_rand_forest_pred(clf, data, result_path, null_dist=None):
    """Makes gene predictions using a previously trained random forest.

    Parameters
    ----------
    clf : GenericClassifier
        random forest sub-class of GenericClassifier
    data : pd.DataFrame
        data frame containing feature information
    result_path : str
        path to save text file result

    Returns
    -------
    tmp_df : pd.DataFrame
        random forest results (already saved to file)
    """
    # perform prediction
    onco_prob, tsg_prob, other_prob = clf.predict()
    true_class = clf.y

    # save features/prediction results
    tmp_df = data.copy()
    gene_order = clf.y.index  # get the order of genes predicted on
    tmp_df['oncogene score'] = pd.Series(onco_prob, index=gene_order)
    tmp_df['tsg score'] = pd.Series(tsg_prob, index=gene_order)
    tmp_df['other score'] = pd.Series(other_prob, index=gene_order)
    tmp_df['driver score'] = 1 - tmp_df['other score']
    pred_class = tmp_df[['other score', 'oncogene score', 'tsg score']].values.argmax(axis=1)
    tmp_df['majority vote class'] = pred_class
    tmp_df['majority vote cancer gene'] = (tmp_df['driver score'] > .5).astype(int)
    tmp_df['training list class'] = true_class
    tmp_df = tmp_df.fillna(0)
    tmp_df = tmp_df.sort(['driver score',], ascending=False)

    if null_dist is not None:
        # add oncogene p-value
        onco_score = tmp_df['oncogene score'].copy()
        onco_score.sort(ascending=False)
        tmp_df['oncogene p-value'] = compute_p_value(onco_score,
                                                     null_dist['oncogene p-value'].dropna())
        tmp_df['oncogene q-value'] = _utils.bh_fdr(tmp_df['oncogene p-value'])

        # add tsg p-value
        tsg_score = tmp_df['tsg score'].copy()
        tsg_score.sort(ascending=False)
        tmp_df['tsg p-value'] = compute_p_value(tsg_score,
                                                null_dist['tsg p-value'].dropna())
        tmp_df['tsg q-value'] = _utils.bh_fdr(tmp_df['tsg p-value'])

        # add driver p-value
        driver_score = tmp_df['driver score'].copy()
        tmp_df['driver p-value'] = compute_p_value(driver_score,
                                                   null_dist['driver p-value'].dropna())
        tmp_df['driver q-value'] = _utils.bh_fdr(tmp_df['driver p-value'])

    tmp_df.to_csv(result_path, sep='\t')

    return tmp_df


def main(cli_opts):
    cfg_opts = _utils.get_output_config('classifier')
    in_opts = _utils.get_input_config('classifier')
    minimum_ct = cli_opts['min_count']

    # get path to features used for classification
    if cli_opts['features']:
        feature_path = cli_opts['features']
    else:
        feature_path = _utils.save_dir + in_opts['gene_feature']

    # read in null distribution p-values
    if not cli_opts['simulated'] and cli_opts['null_distribution']:
        null_pvals = pd.read_csv(cli_opts['null_distribution'], sep='\t',
                                 index_col=0)
    else:
        null_pvals = None

    # use trained classifier if provided
    if cli_opts['trained_classifier']:
        # read in features
        df = pd.read_csv(feature_path, sep='\t', index_col=0)

        logger.info('Running R\'s Random forest . . .')

        # initialize R's random forest
        rrclf = RRandomForest(df,
                              other_sample_ratio=cli_opts['other_ratio'],
                              driver_sample=cli_opts['driver_rate'],
                              ntrees=cli_opts['ntrees'],
                              min_ct=minimum_ct)
        rrclf.clf.load(cli_opts['trained_classifier'])

        # do classification
        pred_results_path = _utils.clf_result_dir + cfg_opts['rrand_forest_pred']
        logger.info('Saving results to {0}'.format(pred_results_path))

        result_df = trained_rand_forest_pred(rrclf, df, pred_results_path, null_pvals)

        if cli_opts['simulated']:
            # driver scores
            driver_score_cts = result_df['driver score'].value_counts()
            driver_score_cts = driver_score_cts.sort_index(ascending=False)
            driver_score_cum_cts = driver_score_cts.cumsum()
            driver_score_pvals = driver_score_cum_cts / float(driver_score_cts.sum())

            # oncogene scores
            onco_score_cts = result_df['oncogene score'].value_counts()
            onco_score_cts = onco_score_cts.sort_index(ascending=False)
            onco_score_cum_cts = onco_score_cts.cumsum()
            onco_score_pvals = onco_score_cum_cts / float(onco_score_cts.sum())

            # tsg score
            tsg_score_cts = result_df['tsg score'].value_counts()
            tsg_score_cts = tsg_score_cts.sort_index(ascending=False)
            tsg_score_cum_cts = tsg_score_cts.cumsum()
            tsg_score_pvals = tsg_score_cum_cts / float(tsg_score_cts.sum())

            # construct null p-value score distribution
            score_ix = set(driver_score_pvals.index) | set(onco_score_pvals.index) | set(tsg_score_pvals.index)
            score_pvals = pd.DataFrame(index=list(score_ix))
            score_pvals['oncogene p-value'] = onco_score_pvals
            score_pvals['tsg p-value'] = tsg_score_pvals
            score_pvals['driver p-value'] = driver_score_pvals
            score_pvals = score_pvals.sort_index(ascending=False)

            score_pvals.to_csv(cli_opts['null_distribution'], sep='\t',
                               index_label='score')

        logger.info('Finished classification.')
        return

    # get oncogene info
    #onco_count_df, onco_pct_df, tsg_ct, tsg_pct = num_pred_by_pct_threshold(minimum_ct)
    #onco_count_df = onco_count_df.sort_index()  # sort df by recurrent mutation cts
    #onco_pct_df = onco_pct_df.sort_index()  # sort df by recurrent mutation cts

    # save results
    #onco_count_df.to_csv(_utils.clf_result_dir + cfg_opts['oncogene_parameters_ct'], sep='\t')
    #onco_pct_df.to_csv(_utils.clf_result_dir + cfg_opts['oncogene_parameters_pct'], sep='\t')

    # get the "normal" results from the 20/20 rule, based on
    # gene_feature_matrix.txt (aka last settings for data_analysis command)
    #logger.info('Generating 20/20 rule predictions . . .')
    #result_df = generate_2020_result(.2, .2, minimum_ct)
    # save result
    #result_df.to_csv(_utils.clf_result_dir + cfg_opts['vogelstein_predictions'], sep='\t')
    # plot results
    #result_df['true class'] = result_df['true class'].apply(lambda x: _utils.class_to_label[x])
    #plot_data.prob_kde(result_df, 'oncogene score',
                       #_utils.clf_plot_dir + cfg_opts['onco_score_kde'],
                       #title='Distribution of Oncogene Scores',
                       #xlabel='Oncogene Score')
    #plot_data.prob_kde(result_df, 'tsg score',
                       #_utils.clf_plot_dir + cfg_opts['tsg_score_kde'],
                       #title='Distribution of TSG Scores',
                       #xlabel='TSG Score')
    #logger.info('Finished generating 20/20 rule predictions.')

    # plot results
    #logger.info('Plotting results of 20/20 rule . . .')
    # plot number of predicted oncogenes while varying parameters
    #tmp_save_path = _utils.clf_plot_dir + cfg_opts['number_oncogenes_plot']
    #tmp_title = r"Landscape 2013 Classifier Predicted Oncogenes"
    #tmp_ylabel = 'Number of Oncogenes'
    #tmp_xlabel = 'Number of Mutations Required for Recurrency'
    #plot_data.onco_mutations_parameter(onco_count_df,
                                       #tmp_save_path,
                                       #title=tmp_title,
                                       #ylabel=tmp_ylabel,
                                       #xlabel=tmp_xlabel)
    # plot percentage of vogelstein's oncogenes recovered
    #tmp_title = 'Percentage of Landscape 2013 Oncogenes Recovered'
    #tmp_ylabel = 'Oncogene Recall'
    #tmp_xlabel = 'Number of Mutations Required for Recurrency'
    #tmp_save_path = _utils.clf_plot_dir + cfg_opts['pct_oncogenes_plot']
    #plot_data.onco_mutations_parameter(onco_pct_df,
                                       #tmp_save_path,
                                       #title=tmp_title,
                                       #ylabel=tmp_ylabel,
                                       #xlabel=tmp_xlabel)
    # plot tsg of number of tsg's predicted
    #tmp_title = 'Landscape 2013 Classifier Predicted TSGs'
    #tmp_ylabel = 'Number of Predicted TSGs'
    #tmp_xlabel = 'Minimum Deleterious Mutations'
    #tmp_save_path = _utils.clf_plot_dir + cfg_opts['number_tsg_plot']
    #plot_data.tsg_mutations_parameter(tsg_ct,
                                      #tmp_save_path,
                                      #title=tmp_title,
                                      #ylabel=tmp_ylabel,
                                      #xlabel=tmp_xlabel)
    # plot recall of tsg while varying tsg score threshold
    #tmp_title = 'Percentage of Landscape 2013 Oncogenes Recovered'
    #tmp_ylabel = 'TSG Recall'
    #tmp_xlabel = 'Minimum Deleterious Mutations'
    #tmp_save_path = _utils.clf_plot_dir + cfg_opts['pct_tsg_plot']
    #plot_data.tsg_mutations_parameter(tsg_pct,
                                      #tmp_save_path,
                                      #title=tmp_title,
                                      #ylabel=tmp_ylabel,
                                      #xlabel=tmp_xlabel)

    df = pd.read_csv(feature_path, sep='\t', index_col=0)

    # plot the 20/20 rule scores
    #plot_data.vogelstein_score_scatter(df.copy(),
                                       #minimum_ct,
                                       #_utils.clf_plot_dir + cfg_opts['2020_score_plot'])
    logger.info('Finished plotting results of 20/20 rule.')

    # R's random forest
    logger.info('Running R\'s Random forest . . .')
    # initialize R's random forest
    rrclf = RRandomForest(df,
                          other_sample_ratio=cli_opts['other_ratio'],
                          driver_sample=cli_opts['driver_rate'],
                          ntrees=cli_opts['ntrees'],
                          min_ct=minimum_ct)

    # analyze classification metrics
    rrclf.kfold_validation()
    rrclf_onco_tpr, rrclf_onco_fpr, rrclf_onco_mean_roc_auc = rrclf.get_onco_roc_metrics()
    rrclf_onco_precision, rrclf_onco_recall, rrclf_onco_mean_pr_auc = rrclf.get_onco_pr_metrics()
    rrclf_tsg_tpr, rrclf_tsg_fpr, rrclf_tsg_mean_roc_auc = rrclf.get_tsg_roc_metrics()
    rrclf_tsg_precision, rrclf_tsg_recall, rrclf_tsg_mean_pr_auc = rrclf.get_tsg_pr_metrics()
    rrclf_driver_precision, rrclf_driver_recall, rrclf_driver_mean_pr_auc = rrclf.get_driver_pr_metrics()
    rrclf_driver_tpr, rrclf_driver_fpr, rrclf_driver_mean_roc_auc = rrclf.get_driver_roc_metrics()

    # plot feature importance
    mean_df = rrclf.mean_importance
    std_df = rrclf.std_importance
    feat_path = _utils.clf_plot_dir + cfg_opts['r_feature_importance_plot']
    plot_data.feature_importance_barplot(mean_df, std_df, feat_path)

    # run predictions using R's random forest
    pred_results_path = _utils.clf_result_dir + cfg_opts['rrand_forest_pred']
    result_df = rand_forest_pred(rrclf, df, result_path=pred_results_path,
                                 null_dist=null_pvals)

    # save a list of oncogenes/tsgs in separate files
    pred_onco = result_df[result_df['predicted class']==_utils.onco_label].index.to_series()
    novel_onco = result_df[(result_df['predicted class']==_utils.onco_label) & (result_df['true class']!=_utils.onco_label)].index.to_series()
    pred_tsg = result_df[result_df['predicted class']==_utils.tsg_label].index.to_series()
    novel_tsg = result_df[(result_df['predicted class']==_utils.tsg_label) & (result_df['true class']!=_utils.tsg_label)].index.to_series()
    pred_driver = result_df[result_df['predicted cancer gene']==1].index.to_series()
    pred_onco.to_csv(_utils.clf_result_dir + cfg_opts['rrf_onco'], sep='\t', index=False, header=None)
    novel_onco.to_csv(_utils.clf_result_dir + cfg_opts['rrf_novel_onco'], sep='\t', index=False, header=None)
    pred_tsg.to_csv(_utils.clf_result_dir + cfg_opts['rrf_tsg'], sep='\t', index=False, header=None)
    novel_tsg.to_csv(_utils.clf_result_dir + cfg_opts['rrf_novel_tsg'], sep='\t', index=False, header=None)
    log_str = ('R Random forest: {0} ({1} novel) oncogenes, '
               '{2} ({3} novel) tsg'.format(len(pred_onco), len(novel_onco),
                                            len(pred_tsg), len(novel_tsg)))
    logger.info(log_str)

    # plot r random forest results
    plot_data.prob_scatter(result_df,
                           plot_path=_utils.clf_plot_dir + cfg_opts['rrand_forest_plot'],
                           title='Sub-sampled Random Forest Predictions')
    plot_data.prob_kde(result_df,
                       col_name='oncogene score',
                       save_path=_utils.clf_plot_dir + cfg_opts['onco_kde_rrand_forest'],
                       title='Distribution of Oncogene Scores (sub-sampled random forest)')
    plot_data.prob_kde(result_df,
                       col_name='tsg score',
                       save_path=_utils.clf_plot_dir + cfg_opts['tsg_kde_rrand_forest'],
                       title='Distribution of TSG Scores (sub-sampled random forest)')
    #result_df['driver gene probability'] = result_df['oncogene probability'] + result_df['tsg probability']
    plot_data.prob_kde(result_df,
                       col_name='driver score',
                       save_path='results/classify/plots/r_random_forest_driver_prob.kde.png',
                       title='Distribution of Driver Gene Score (sub-sampled random forest)')
    #plot_data.sample_boxplot(pred_onco,
                             #pred_tsg,
                             #pred_driver,
                             #save_path_type=_utils.clf_plot_dir + cfg_opts['rrf_sample_pct_type_boxplot'],
                             #save_path_driver=_utils.clf_plot_dir + cfg_opts['rrf_sample_pct_driver_boxplot'])
    logger.info('Finished running sub-sampled Random Forest')

    # scikit learns' random forest
    logger.info('Running Random forest . . .')
    rclf = RandomForest(df, ntrees=cli_opts['ntrees'], min_ct=minimum_ct)
    rclf.kfold_validation()
    rclf_onco_tpr, rclf_onco_fpr, rclf_onco_mean_roc_auc = rclf.get_onco_roc_metrics()
    rclf_onco_precision, rclf_onco_recall, rclf_onco_mean_pr_auc = rclf.get_onco_pr_metrics()
    rclf_tsg_tpr, rclf_tsg_fpr, rclf_tsg_mean_roc_auc = rclf.get_tsg_roc_metrics()
    rclf_tsg_precision, rclf_tsg_recall, rclf_tsg_mean_pr_auc = rclf.get_tsg_pr_metrics()
    rclf_driver_precision, rclf_driver_recall, rclf_driver_mean_pr_auc = rclf.get_driver_pr_metrics()
    rclf_driver_tpr, rclf_driver_fpr, rclf_driver_mean_roc_auc = rclf.get_driver_roc_metrics()

    # plot feature importance
    mean_df = rclf.mean_importance
    std_df = rclf.std_importance
    feat_path = _utils.clf_plot_dir + cfg_opts['feature_importance_plot']
    plot_data.feature_importance_barplot(mean_df, std_df, feat_path)

    # predict using scikit learn's random forest
    pred_path = _utils.clf_result_dir + cfg_opts['rand_forest_pred']
    result_df = rand_forest_pred(rclf, df,
                                 result_path=pred_path,
                                 null_dist=null_pvals)

    # save a list of oncogenes/tsgs in separate files
    pred_onco = result_df[result_df['predicted class']==_utils.onco_label].index.to_series()
    novel_onco = result_df[(result_df['predicted class']==_utils.onco_label) & (result_df['true class']!=_utils.onco_label)].index.to_series()
    pred_tsg = result_df[result_df['predicted class']==_utils.tsg_label].index.to_series()
    novel_tsg = result_df[(result_df['predicted class']==_utils.tsg_label) & (result_df['true class']!=_utils.tsg_label)].index.to_series()
    pred_onco.to_csv(_utils.clf_result_dir + cfg_opts['rf_onco'], sep='\t', index=False, header=None)
    novel_onco.to_csv(_utils.clf_result_dir + cfg_opts['rf_novel_onco'], sep='\t', index=False, header=None)
    pred_tsg.to_csv(_utils.clf_result_dir + cfg_opts['rf_tsg'], sep='\t', index=False, header=None)
    novel_tsg.to_csv(_utils.clf_result_dir + cfg_opts['rf_novel_tsg'], sep='\t', index=False, header=None)
    log_str = ('Random forest: {0} ({1} novel) oncogenes, '
               '{2} ({3} novel) tsg'.format(len(pred_onco), len(novel_onco),
                                            len(pred_tsg), len(novel_tsg)))
    logger.info(log_str)

    # plot random forest result
    plot_data.prob_scatter(result_df,
                           plot_path=_utils.clf_plot_dir + cfg_opts['rand_forest_plot'],
                           title='Random Forest Predictions')
    plot_data.prob_kde(result_df,
                       col_name='oncogene score',
                       save_path=_utils.clf_plot_dir + cfg_opts['onco_kde_rand_forest'],
                       title='Distribution of Oncogene Probabilities (random forest)')
    plot_data.prob_kde(result_df,
                       col_name='tsg score',
                       save_path=_utils.clf_plot_dir + cfg_opts['tsg_kde_rand_forest'],
                       title='Distribution of TSG Score (random forest)')
    # result_df['driver gene probability'] = result_df['oncogene probability'] + result_df['tsg probability']
    plot_data.prob_kde(result_df,
                       col_name='driver score',
                       save_path='results/classify/plots/random_forest_driver_prob.kde.png',
                       title='Distribution of Driver Gene Scores (random forest)')
    logger.info('Finished running Random Forest')

    # multinomial naive bayes
    logger.info('Running Naive Bayes . . .')
    #nbclf = MultinomialNaiveBayes(df, min_ct=minimum_ct)
    #nbclf.kfold_validation()
    #nbclf_onco_tpr, nbclf_onco_fpr, nbclf_onco_mean_roc_auc = nbclf.get_onco_roc_metrics()
    #nbclf_onco_precision, nbclf_onco_recall, nbclf_onco_mean_pr_auc = nbclf.get_onco_pr_metrics()
    #nbclf_tsg_tpr, nbclf_tsg_fpr, nbclf_tsg_mean_roc_auc = nbclf.get_tsg_roc_metrics()
    #nbclf_tsg_precision, nbclf_tsg_recall, nbclf_tsg_mean_pr_auc = nbclf.get_tsg_pr_metrics()
    #nbclf_driver_precision, nbclf_driver_recall, nbclf_driver_mean_pr_auc = nbclf.get_driver_pr_metrics()
    #nbclf_driver_tpr, nbclf_driver_fpr, nbclf_driver_mean_roc_auc = nbclf.get_driver_roc_metrics()
    logger.info('Finished Naive Bayes.')

    # dummy classifier, predict most frequent
    logger.info('Running Dummy Classifier. . .')
    dclf = DummyClf(df,
                    strategy='most_frequent',
                    min_ct=minimum_ct,
                    weight=False)
    dclf.kfold_validation()
    dclf_onco_tpr, dclf_onco_fpr, dclf_onco_mean_roc_auc = dclf.get_onco_roc_metrics()
    dclf_onco_precision, dclf_onco_recall, dclf_onco_mean_pr_auc = dclf.get_onco_pr_metrics()
    dclf_tsg_tpr, dclf_tsg_fpr, dclf_tsg_mean_roc_auc = dclf.get_tsg_roc_metrics()
    dclf_tsg_precision, dclf_tsg_recall, dclf_tsg_mean_pr_auc = dclf.get_tsg_pr_metrics()
    dclf_driver_tpr, dclf_driver_fpr, dclf_driver_mean_roc_auc = dclf.get_driver_roc_metrics()
    logger.info('Finished dummy classifier.')

    # plot oncogene roc figure
    random_forest_str = 'random forest (AUC = %0.3f)' % rclf_onco_mean_roc_auc
    rrandom_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_onco_mean_roc_auc
    #naive_bayes_str = 'naive bayes (AUC = %0.3f)' % nbclf_onco_mean_roc_auc
    dummy_str = 'dummy (AUC = %0.3f)' % dclf_onco_mean_roc_auc
    rclf_onco_mean_tpr = np.mean(rclf_onco_tpr, axis=0)
    rrclf_onco_mean_tpr = np.mean(rrclf_onco_tpr, axis=0)
    #nbclf_onco_mean_tpr = np.mean(nbclf_onco_tpr, axis=0)
    dclf_onco_mean_tpr = np.mean(dclf_onco_tpr, axis=0)
    df = pd.DataFrame({random_forest_str: rclf_onco_mean_tpr,
                       rrandom_forest_str: rrclf_onco_mean_tpr,
                       #naive_bayes_str: nbclf_onco_mean_tpr,
                       dummy_str: dclf_onco_mean_tpr},
                      index=rclf_onco_fpr)
    line_style = {dummy_str: '--',
                  random_forest_str: '-',
                  rrandom_forest_str: '-',
                  #naive_bayes_str: '-'
                  }
    save_path = _utils.clf_plot_dir + cfg_opts['roc_plot_oncogene']
    plot_data.receiver_operator_curve(df, save_path, line_style)

    # plot tsg roc figure
    random_forest_str = 'random forest (AUC = %0.3f)' % rclf_tsg_mean_roc_auc
    r_random_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_tsg_mean_roc_auc
    #naive_bayes_str = 'naive bayes (AUC = %0.3f)' % nbclf_tsg_mean_roc_auc
    dummy_str = 'dummy (AUC = %0.3f)' % dclf_tsg_mean_roc_auc
    rclf_tsg_mean_tpr = np.mean(rclf_tsg_tpr, axis=0)
    rrclf_tsg_mean_tpr = np.mean(rrclf_tsg_tpr, axis=0)
    #nbclf_tsg_mean_tpr = np.mean(nbclf_tsg_tpr, axis=0)
    dclf_tsg_mean_tpr = np.mean(dclf_tsg_tpr, axis=0)
    df = pd.DataFrame({random_forest_str: rclf_tsg_mean_tpr,
                       r_random_forest_str: rrclf_tsg_mean_tpr,
                       #naive_bayes_str: nbclf_tsg_mean_tpr,
                       dummy_str: dclf_tsg_mean_tpr},
                      index=rclf_tsg_fpr)
    line_style = {dummy_str: '--',
                  random_forest_str: '-',
                  #naive_bayes_str: '-',
                  }
    save_path = _utils.clf_plot_dir + cfg_opts['roc_plot_tsg']
    plot_data.receiver_operator_curve(df, save_path, line_style)

    # plot driver roc figure
    random_forest_str = 'random forest (AUC = %0.3f)' % rclf_driver_mean_roc_auc
    r_random_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_driver_mean_roc_auc
    #naive_bayes_str = 'naive bayes (AUC = %0.3f)' % nbclf_driver_mean_roc_auc
    dummy_str = 'dummy (AUC = %0.3f)' % dclf_driver_mean_roc_auc
    rclf_driver_mean_tpr = np.mean(rclf_driver_tpr, axis=0)
    rrclf_driver_mean_tpr = np.mean(rrclf_driver_tpr, axis=0)
    #nbclf_driver_mean_tpr = np.mean(nbclf_driver_tpr, axis=0)
    dclf_driver_mean_tpr = np.mean(dclf_driver_tpr, axis=0)
    df = pd.DataFrame({random_forest_str: rclf_driver_mean_tpr,
                       r_random_forest_str: rrclf_driver_mean_tpr,
                       #naive_bayes_str: nbclf_driver_mean_tpr,
                       dummy_str: dclf_driver_mean_tpr},
                      index=rclf_driver_fpr)
    line_style = {dummy_str: '--',
                  random_forest_str: '-',
                  #naive_bayes_str:'-',
                  }
    save_path = _utils.clf_plot_dir + cfg_opts['roc_plot_driver']
    plot_data.receiver_operator_curve(df, save_path, line_style)

    # plot oncogene pr figure
    random_forest_str = 'random forest (AUC = %0.3f)' % rclf_onco_mean_pr_auc
    rrandom_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_onco_mean_pr_auc
    #naive_bayes_str = 'naive bayes (AUC = %0.3f)' % nbclf_onco_mean_pr_auc
    dummy_str = 'dummy (AUC = %0.3f)' % dclf_onco_mean_pr_auc
    rclf_onco_mean_precision = np.mean(rclf_onco_precision, axis=0)
    rrclf_onco_mean_precision = np.mean(rrclf_onco_precision, axis=0)
    #nbclf_onco_mean_precision = np.mean(nbclf_onco_precision, axis=0)
    dclf_onco_mean_precision = np.mean(dclf_onco_precision, axis=0)
    df = pd.DataFrame({random_forest_str: rclf_onco_mean_precision,
                       rrandom_forest_str: rrclf_onco_mean_precision,
                       #naive_bayes_str: nbclf_onco_mean_precision
                       },
                      index=rclf_onco_recall)
    #rclf_onco_sem_precision = stats.sem(rclf_onco_precision, axis=0)
    #nbclf_onco_sem_precision = stats.sem(nbclf_onco_precision, axis=0)
    #dclf_onco_sem_precision = stats.sem(dclf_onco_precision, axis=0)
    #sem_df = df.copy()
    #sem_df[random_forest_str] = rclf_onco_sem_precision
    #sem_df[naive_bayes_str] = nbclf_onco_sem_precision
    #sem_df[dummy_str] = dclf_onco_sem_precision
    line_style = {dummy_str: '--',
                  random_forest_str: '-',
                  rrandom_forest_str: '-',
                  #naive_bayes_str:'-'
                  }
    save_path = _utils.clf_plot_dir + cfg_opts['pr_plot_oncogene']
    plot_data.precision_recall_curve(df, save_path, line_style,
                                     #sem_df,
                                     title='Oncogene Precision-Recall Curve')

    # plot tsg pr figure
    random_forest_str = 'random forest (AUC = %0.3f)' % rclf_tsg_mean_pr_auc
    r_random_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_tsg_mean_pr_auc
    #naive_bayes_str = 'naive bayes (AUC = %0.3f)' % nbclf_tsg_mean_pr_auc
    dummy_str = 'dummy (AUC = %0.3f)' % dclf_tsg_mean_pr_auc
    rclf_tsg_mean_precision = np.mean(rclf_tsg_precision, axis=0)
    rrclf_tsg_mean_precision = np.mean(rrclf_tsg_precision, axis=0)
    #nbclf_tsg_mean_precision = np.mean(nbclf_tsg_precision, axis=0)
    dclf_tsg_mean_precision = np.mean(dclf_tsg_precision, axis=0)
    df = pd.DataFrame({random_forest_str: rclf_tsg_mean_precision,
                       r_random_forest_str: rrclf_tsg_mean_precision,
                       #naive_bayes_str: nbclf_tsg_mean_precision
                       },
                      index=rclf_tsg_recall)
    line_style = {dummy_str: '--',
                  random_forest_str: '-',
                  r_random_forest_str: '-',
                  #naive_bayes_str:'-'
                  }
    save_path = _utils.clf_plot_dir + cfg_opts['pr_plot_tsg']
    plot_data.precision_recall_curve(df, save_path, line_style,
                                     title='TSG Precision-Recall Curve')


    # plot driver gene pr figure
    random_forest_str = 'random forest (AUC = %0.3f)' % rclf_driver_mean_pr_auc
    r_random_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_driver_mean_pr_auc
    # naive_bayes_str = 'naive bayes (AUC = %0.3f)' % nbclf_driver_mean_pr_auc
    rclf_driver_mean_precision = np.mean(rclf_driver_precision, axis=0)
    rrclf_driver_mean_precision = np.mean(rrclf_driver_precision, axis=0)
    #nbclf_driver_mean_precision = np.mean(nbclf_driver_precision, axis=0)
    df = pd.DataFrame({random_forest_str: rclf_driver_mean_precision,
                       r_random_forest_str: rrclf_driver_mean_precision,
                       # naive_bayes_str: nbclf_driver_mean_precision
                       },
                      index=rclf_driver_recall)
    line_style = {dummy_str: '--',
                  random_forest_str: '-',
                  r_random_forest_str: '-',
                  #naive_bayes_str:'-'
                  }
    save_path = _utils.clf_plot_dir + cfg_opts['pr_plot_driver']
    plot_data.precision_recall_curve(df, save_path, line_style,
                                     title='Driver Precision-Recall Curve')


if __name__ == "__main__":
    main()
