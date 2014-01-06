from __future__ import division
from vogelstein_classifier import VogelsteinClassifier
from random_forest_clf import RandomForest
from multinomial_nb_clf import MultinomialNaiveBayes
from dummy_clf import DummyClf
from r_random_forest_clf import RRandomForest
import utils.python.util as _utils
import plot_data
import pandas as pd
import numpy as np
import glob
import re
import utils.python.plot as myplt
import logging

logger = logging.getLogger(__name__)

def calc_onco_info(df, onco_pct, tsg_pct, min_ct):
    # calculate the number of genes classified as oncogene
    vclf = VogelsteinClassifier(onco_pct, tsg_pct, min_count=min_ct)
    df['total'] = df.T.sum()
    input_list = ((row['recurrent missense'] + row['recurrent indel'],
                   row['frame shift'] + row['nonsense'] + row['lost stop'] + row['no protein'],
                   row['total'])
                  for i, row in df.iterrows())
    df['2020 class'] = vclf.predict_list(input_list)
    class_cts = df['2020 class'].value_counts()

    # calculate the pct of known oncogenes found
    df['curated class'] = [_utils.classify_gene(gene)
                           for gene in df.index.tolist()]
    tmpdf = df.copy()  # prevent wierd behavior
    known_onco = tmpdf[tmpdf['curated class']=='oncogene']
    num_onco_found = len(known_onco[known_onco['2020 class']=='oncogene'])
    total_onco = len(known_onco)  # total number of oncogenes with counts
    pct_onco_found = num_onco_found / total_onco

    return class_cts['oncogene'], pct_onco_found


def num_onco_by_recurrent_mutations(onco_pct, tsg_pct, min_ct):
    """Count number of oncogenes while varying the definition of recurrency"""
    # calculate counts for oncogenes/tsg with varying the required the number
    # of mutations to define a recurrent position
    file_match_pattern = './data_analysis/results/genes/gene_feature_matrix.r*.txt'
    gene_design_matrix_paths = glob.glob(file_match_pattern)
    onco_ct_list, onco_pct_list = [], []  # list of cts/pct for oncogenes
    for file_path in gene_design_matrix_paths:
        tmp_df = pd.read_csv(file_path, sep='\t', index_col=0)
        tmp_ct, tmp_pct = calc_onco_info(tmp_df,
                                         onco_pct=onco_pct,  # pct thresh for onco
                                         tsg_pct=tsg_pct,  # pct thresh for tsg
                                         min_ct=min_ct)  # min count for a gene
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


def num_onco_by_pct_threshold(min_ct):
    # initialization of dataframe
    cts, pct = num_onco_by_recurrent_mutations(.2, .2, min_ct)
    df_ct = pd.DataFrame(index=cts.index)
    df_pct = pd.DataFrame(index=pct.index)

    # test different percentage thresholds
    for threshold in np.arange(.15, .5, .05):
        tmp_ct, tmp_pct = num_onco_by_recurrent_mutations(threshold, threshold, min_ct)
        df_ct[str(threshold)] = tmp_ct
        df_pct[str(threshold)] = tmp_pct
    return df_ct, df_pct


def rand_forest_pred(clf, data, result_path,
                     plot_path, plot_title):
    """Makes gene predictions using a random forest classifier."""
    onco_prob, tsg_prob, other_prob = clf.kfold_prediction()
    true_class = clf.y
    tmp_df = data.copy()
    tmp_df['true class'] = true_class
    # tmp_df['predicted class'] = pred
    tmp_df['onco prob class'] = onco_prob
    tmp_df['tsg prob class'] = tsg_prob
    tmp_df['other prob class'] = other_prob
    pred_class = tmp_df[['other prob class', 'onco prob class', 'tsg prob class']].values.argmax(axis=1)
    tmp_df['predicted class'] = pred_class
    tmp_df = tmp_df.fillna(0)
    tmp_df = tmp_df.sort(['onco prob class', 'tsg prob class'], ascending=False)
    tmp_df.to_csv(result_path, sep='\t')
    myplt.scatter(tmp_df['onco prob class'],
                  tmp_df['tsg prob class'],
                  plot_path,
                  xlabel='Oncogene Probability',
                  ylabel='TSG Probability',
                  title='R\'s Random Forest Predictions',
                  colors='#348ABD')


def main(cli_opts):
    cfg_opts = _utils.get_output_config('classifier')
    in_opts = _utils.get_input_config('classifier')
    minimum_ct = cli_opts['min_count']

    # get oncogene info
    count_df, pct_df = num_onco_by_pct_threshold(minimum_ct)
    count_df = count_df.sort_index()  # sort df by recurrent mutation cts
    pct_df = pct_df.sort_index()  # sort df by recurrent mutation cts

    # save results
    count_df.to_csv(_utils.clf_result_dir + cfg_opts['oncogene_parameters_ct'], sep='\t')
    pct_df.to_csv(_utils.clf_result_dir + cfg_opts['oncogene_parameters_pct'], sep='\t')

    # plot results
    logger.info('Plotting results of 20/20 rule . . .')
    # plot number of predicted oncogenes while varying parameters
    tmp_save_path = _utils.clf_plot_dir + cfg_opts['number_oncogenes_plot']
    tmp_title = r"Vogelstein's Classifier Predicted Oncogenes"
    tmp_ylabel = 'Number of Oncogenes'
    tmp_xlabel = 'Number of Mutations Required for Recurrency'
    plot_data.onco_mutations_parameter(count_df,
                                       tmp_save_path,
                                       title=tmp_title,
                                       ylabel=tmp_ylabel,
                                       xlabel=tmp_xlabel)
    # plot percentage of vogelstein's oncogenes recovered
    tmp_title = 'Percentage of Vogelstein\'s Oncogenes Recovered'
    tmp_ylabel = 'Oncogene Recall'
    tmp_xlabel = 'Number of Mutations Required for Recurrency'
    tmp_save_path = _utils.clf_plot_dir + cfg_opts['pct_oncogenes_plot']
    plot_data.onco_mutations_parameter(pct_df,
                                       tmp_save_path,
                                       title=tmp_title,
                                       ylabel=tmp_ylabel,
                                       xlabel=tmp_xlabel)

    df = pd.read_csv(in_opts['gene_feature'],
                     sep='\t', index_col=0)

    # plot the 20/20 rule scores
    plot_data.vogelstein_score_scatter(df.copy(),
                                       minimum_ct,
                                       _utils.clf_plot_dir + cfg_opts['2020_score_plot'])
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

    # run predictions using R's random forest
    rand_forest_pred(rrclf, df,
                     result_path=_utils.clf_result_dir + cfg_opts['rrand_forest_pred'],
                     plot_path=_utils.clf_plot_dir + cfg_opts['rrand_forest_plot'],
                     plot_title='R\'s Random Forest Predictions')
    logger.info('Finished running R\'s Random Forest')

    # scikit learns' random forest
    logger.info('Running Random forest . . .')
    rclf = RandomForest(df, ntrees=cli_opts['ntrees'], min_ct=minimum_ct)
    rclf.kfold_validation()
    rclf_onco_tpr, rclf_onco_fpr, rclf_onco_mean_roc_auc = rclf.get_onco_roc_metrics()
    rclf_onco_precision, rclf_onco_recall, rclf_onco_mean_pr_auc = rclf.get_onco_pr_metrics()
    rclf_tsg_tpr, rclf_tsg_fpr, rclf_tsg_mean_roc_auc = rclf.get_tsg_roc_metrics()
    rclf_tsg_precision, rclf_tsg_recall, rclf_tsg_mean_pr_auc = rclf.get_tsg_pr_metrics()

    # plot feature importance
    mean_df = rclf.mean_importance
    std_df = rclf.std_importance
    plot_data.feature_importance_barplot(mean_df,
                                         std_df,
                                         _utils.clf_plot_dir + cfg_opts['feature_importance_plot'])

    # predict using scikit learn's random forest
    rand_forest_pred(rclf, df,
                     result_path=_utils.clf_result_dir + cfg_opts['rand_forest_pred'],
                     plot_path=_utils.clf_plot_dir + cfg_opts['rand_forest_plot'],
                     plot_title='Random Forest Predictions')
    logger.info('Finished running Random Forest')

    # multinomial naive bayes
    logger.info('Running Naive Bayes . . .')
    nbclf = MultinomialNaiveBayes(df, min_ct=minimum_ct)
    nbclf.kfold_validation()
    nbclf_onco_tpr, nbclf_onco_fpr, nbclf_onco_mean_roc_auc = nbclf.get_onco_roc_metrics()
    nbclf_onco_precision, nbclf_onco_recall, nbclf_onco_mean_pr_auc = nbclf.get_onco_pr_metrics()
    nbclf_tsg_tpr, nbclf_tsg_fpr, nbclf_tsg_mean_roc_auc = nbclf.get_tsg_roc_metrics()
    nbclf_tsg_precision, nbclf_tsg_recall, nbclf_tsg_mean_pr_auc = nbclf.get_tsg_pr_metrics()
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
    logger.info('Finished dummy classifier.')

    # plot oncogene roc figure
    random_forest_str = 'random forest (AUC = %0.3f)' % rclf_onco_mean_roc_auc
    rrandom_forest_str = 'R\'s random forest (AUC = %0.3f)' % rrclf_onco_mean_roc_auc
    naive_bayes_str = 'naive bayes (AUC = %0.3f)' % nbclf_onco_mean_roc_auc
    dummy_str = 'dummy (AUC = %0.3f)' % dclf_onco_mean_roc_auc
    rclf_onco_mean_tpr = np.mean(rclf_onco_tpr, axis=0)
    rrclf_onco_mean_tpr = np.mean(rrclf_onco_tpr, axis=0)
    nbclf_onco_mean_tpr = np.mean(nbclf_onco_tpr, axis=0)
    dclf_onco_mean_tpr = np.mean(dclf_onco_tpr, axis=0)
    df = pd.DataFrame({random_forest_str: rclf_onco_mean_tpr,
                       naive_bayes_str: nbclf_onco_mean_tpr,
                       rrandom_forest_str: rrclf_onco_mean_tpr,
                       dummy_str: dclf_onco_mean_tpr},
                      index=rclf_onco_fpr)
    line_style = {dummy_str: '--',
                  random_forest_str: '-',
                  rrandom_forest_str: '-',
                  naive_bayes_str:'-'}
    save_path = _utils.clf_plot_dir + cfg_opts['roc_plot_oncogene']
    plot_data.receiver_operator_curve(df, save_path, line_style)

    # plot tsg roc figure
    random_forest_str = 'random forest (AUC = %0.3f)' % rclf_tsg_mean_roc_auc
    naive_bayes_str = 'naive bayes (AUC = %0.3f)' % nbclf_tsg_mean_roc_auc
    dummy_str = 'dummy (AUC = %0.3f)' % dclf_tsg_mean_roc_auc
    rclf_tsg_mean_tpr = np.mean(rclf_tsg_tpr, axis=0)
    nbclf_tsg_mean_tpr = np.mean(nbclf_tsg_tpr, axis=0)
    dclf_tsg_mean_tpr = np.mean(dclf_tsg_tpr, axis=0)
    df = pd.DataFrame({random_forest_str: rclf_tsg_mean_tpr,
                       naive_bayes_str: nbclf_tsg_mean_tpr,
                       dummy_str: dclf_tsg_mean_tpr},
                      index=rclf_tsg_fpr)
    line_style = {dummy_str: '--',
                  random_forest_str: '-',
                  naive_bayes_str:'-'}
    save_path = _utils.clf_plot_dir + cfg_opts['roc_plot_tsg']
    plot_data.receiver_operator_curve(df, save_path, line_style)

    # plot oncogene pr figure
    random_forest_str = 'random forest (AUC = %0.3f)' % rclf_onco_mean_pr_auc
    rrandom_forest_str = 'R\'s random forest (AUC = %0.3f)' % rrclf_onco_mean_pr_auc
    naive_bayes_str = 'naive bayes (AUC = %0.3f)' % nbclf_onco_mean_pr_auc
    dummy_str = 'dummy (AUC = %0.3f)' % dclf_onco_mean_pr_auc
    rclf_onco_mean_precision = np.mean(rclf_onco_precision, axis=0)
    rrclf_onco_mean_precision = np.mean(rrclf_onco_precision, axis=0)
    nbclf_onco_mean_precision = np.mean(nbclf_onco_precision, axis=0)
    dclf_onco_mean_precision = np.mean(dclf_onco_precision, axis=0)
    df = pd.DataFrame({random_forest_str: rclf_onco_mean_precision,
                       rrandom_forest_str: rrclf_onco_mean_precision,
                       naive_bayes_str: nbclf_onco_mean_precision,
                       dummy_str: dclf_onco_mean_precision},
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
                  naive_bayes_str:'-'}
    save_path = _utils.clf_plot_dir + cfg_opts['pr_plot_oncogene']
    plot_data.precision_recall_curve(df, save_path, line_style,
                                     #sem_df,
                                     title='Oncogene Precision-Recall Curve')

    # plot tsg pr figure
    random_forest_str = 'random forest (AUC = %0.3f)' % rclf_tsg_mean_pr_auc
    naive_bayes_str = 'naive bayes (AUC = %0.3f)' % nbclf_tsg_mean_pr_auc
    dummy_str = 'dummy (AUC = %0.3f)' % dclf_tsg_mean_pr_auc
    rclf_tsg_mean_precision = np.mean(rclf_tsg_precision, axis=0)
    nbclf_tsg_mean_precision = np.mean(nbclf_tsg_precision, axis=0)
    dclf_tsg_mean_precision = np.mean(dclf_tsg_precision, axis=0)
    df = pd.DataFrame({random_forest_str: rclf_tsg_mean_precision,
                       naive_bayes_str: nbclf_tsg_mean_precision,
                       dummy_str: dclf_tsg_mean_precision},
                      index=rclf_tsg_recall)
    line_style = {dummy_str: '--',
                  random_forest_str: '-',
                  naive_bayes_str:'-'}
    save_path = _utils.clf_plot_dir + cfg_opts['pr_plot_tsg']
    plot_data.precision_recall_curve(df, save_path, line_style,
                                     title='TSG Precision-Recall Curve')


if __name__ == "__main__":
    main()
