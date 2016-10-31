from __future__ import division
from src.classify.python.vogelstein_classifier import VogelsteinClassifier
from src.classify.python.multinomial_nb_clf import MultinomialNaiveBayes
from src.classify.python.dummy_clf import DummyClf
from src.classify.python.r_random_forest_clf import RRandomForest
import src.utils.python.util as _utils
import pandas as pd
import numpy as np
import glob
import re
import bisect
import logging

# skip plotting if they don't have matplotlib
try:
    import src.classify.python.plot_data as plot_data
except ImportError:
    pass

logger = logging.getLogger(__name__)

def compute_p_value(scores, null_p_values):
    """Get the p-value for each score by examining the list null distribution
    where scores are obtained by a certain probability.
    
    NOTE: uses score2pval function
    
    Parameters
    ----------
    scores : pd.Series
        series of observed scores
    null_p_values: pd.Series
        Empirical null distribution, index are scores and values are p values
    
    Returns
    -------
    pvals : pd.Series
        Series of p values for scores
    """
    num_scores = len(scores)
    pvals = pd.Series(np.zeros(num_scores))
    null_p_val_scores = list(reversed(null_p_values.index.tolist()))
    #null_p_values = null_p_values.ix[null_p_val_scores].copy()
    null_p_values.sort_values(inplace=True, ascending=False)
    pvals = scores.apply(lambda x: score2pval(x, null_p_val_scores, null_p_values))
    return pvals


def score2pval(score, null_scores, null_pvals):
    """Looks up the P value from the empirical null distribution based on the provided
    score.
    
    NOTE: null_scores and null_pvals should be sorted in ascending order.
    
    Parameters
    ----------
    score : float
        score to look up P value for
    null_scores : list
        list of scores that have a non-NA value
    null_pvals : pd.Series
        a series object with the P value for the scores found in null_scores
    
    Returns
    -------
    pval : float
        P value for requested score
    """
    # find position in simulated null distribution
    pos = bisect.bisect_right(null_scores, score)

    # if the score is beyond any simulated values, then report
    # a p-value of zero
    if pos == null_pvals.size and score > null_scores[-1]:
        return 0
    # condition needed to prevent an error
    # simply get last value, if it equals the last value
    elif pos == null_pvals.size:
        return null_pvals.iloc[pos-1]
    # normal case, just report the corresponding p-val from simulations
    else:
        return null_pvals.iloc[pos]


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
    tmp_df = tmp_df.sort_values(by=['driver score',], ascending=False)

    if null_dist is not None:
        # add oncogene p-value
        onco_score = tmp_df['oncogene score'].copy()
        onco_score.sort_values(inplace=True, ascending=False)
        tmp_df['oncogene p-value'] = compute_p_value(onco_score,
                                                     null_dist['oncogene p-value'].dropna())
        tmp_df['oncogene q-value'] = _utils.bh_fdr(tmp_df['oncogene p-value'])

        # add tsg p-value
        tsg_score = tmp_df['tsg score'].copy()
        tsg_score.sort_values(inplace=True, ascending=False)
        tmp_df['tsg p-value'] = compute_p_value(tsg_score,
                                                null_dist['tsg p-value'].dropna())
        tmp_df['tsg q-value'] = _utils.bh_fdr(tmp_df['tsg p-value'])

        # add driver p-values
        driver_score = tmp_df['driver score'].copy()
        driver_score.sort_values(inplace=True, ascending=False)
        tmp_df['driver p-value'] = compute_p_value(driver_score,
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
    null_dist : pd.DataFrame (default: None)
        dataframe relating scores to p values

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
    tmp_df = tmp_df.sort_values(by=['driver score',], ascending=False)

    if null_dist is not None:
        # add oncogene p-value
        onco_score = tmp_df['oncogene score'].copy()
        onco_score.sort_values(inplace=True, ascending=False)
        tmp_df['oncogene p-value'] = compute_p_value(onco_score,
                                                     null_dist['oncogene p-value'].dropna())
        tmp_df['oncogene q-value'] = _utils.bh_fdr(tmp_df['oncogene p-value'])

        # add tsg p-value
        tsg_score = tmp_df['tsg score'].copy()
        tsg_score.sort_values(inplace=True, ascending=False)
        tmp_df['tsg p-value'] = compute_p_value(tsg_score,
                                                null_dist['tsg p-value'].dropna())
        tmp_df['tsg q-value'] = _utils.bh_fdr(tmp_df['tsg p-value'])

        # add driver p-value
        driver_score = tmp_df['driver score'].copy()
        driver_score.sort_values(inplace=True, ascending=False)
        tmp_df['driver p-value'] = compute_p_value(driver_score,
                                                   null_dist['driver p-value'].dropna())
        tmp_df['driver q-value'] = _utils.bh_fdr(tmp_df['driver p-value'])
    else:
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

        logger.info('Running Random forest . . .')

        # initialize R's random forest
        rrclf = RRandomForest(df,
                              other_sample_ratio=cli_opts['other_ratio'],
                              driver_sample=cli_opts['driver_rate'],
                              ntrees=cli_opts['ntrees'],
                              seed=cli_opts['random_seed'])
        rrclf.clf.load(cli_opts['trained_classifier'])

        if cli_opts['simulated']:
            # do classification
            result_df = trained_rand_forest_pred(rrclf, df, None, null_pvals)

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
        else:
            # do classification
            pred_results_path = _utils.clf_result_dir + cfg_opts['rrand_forest_pred']
            logger.info('Saving results to {0}'.format(pred_results_path))
            result_df = trained_rand_forest_pred(rrclf, df, pred_results_path, null_pvals)
            result_df.to_csv(pred_results_path, sep='\t')

        logger.info('Finished classification.')
        return

    df = pd.read_csv(feature_path, sep='\t', index_col=0)

    # R's random forest
    logger.info('Running Random forest . . .')
    # initialize R's random forest
    rrclf = RRandomForest(df,
                          other_sample_ratio=cli_opts['other_ratio'],
                          driver_sample=cli_opts['driver_rate'],
                          ntrees=cli_opts['ntrees'],
                          seed=cli_opts['random_seed'])

    # analyze classification metrics
    rrclf.kfold_validation()
    rrclf_onco_tpr, rrclf_onco_fpr, rrclf_onco_mean_roc_auc = rrclf.get_onco_roc_metrics()
    rrclf_onco_precision, rrclf_onco_recall, rrclf_onco_mean_pr_auc = rrclf.get_onco_pr_metrics()
    rrclf_tsg_tpr, rrclf_tsg_fpr, rrclf_tsg_mean_roc_auc = rrclf.get_tsg_roc_metrics()
    rrclf_tsg_precision, rrclf_tsg_recall, rrclf_tsg_mean_pr_auc = rrclf.get_tsg_pr_metrics()
    rrclf_driver_precision, rrclf_driver_recall, rrclf_driver_mean_pr_auc = rrclf.get_driver_pr_metrics()
    rrclf_driver_tpr, rrclf_driver_fpr, rrclf_driver_mean_roc_auc = rrclf.get_driver_roc_metrics()

    # skip if no matplotlib
    try:
        # plot feature importance
        mean_df = rrclf.mean_importance
        std_df = rrclf.std_importance
        feat_path = _utils.clf_plot_dir + cfg_opts['r_feature_importance_plot']
        plot_data.feature_importance_barplot(mean_df, std_df, feat_path)
    except:
        pass

    # run predictions using R's random forest
    pred_results_path = _utils.clf_result_dir + cfg_opts['rrand_forest_pred']
    result_df = rand_forest_pred(rrclf, df, result_path=pred_results_path,
                                 null_dist=null_pvals)

    # save a list of oncogenes/tsgs in separate files
    if null_pvals is None:
        pred_onco = result_df[result_df['majority vote class']==_utils.onco_label].index.to_series()
        novel_onco = result_df[(result_df['majority vote class']==_utils.onco_label) & (result_df['training list class']!=_utils.onco_label)].index.to_series()
        pred_tsg = result_df[result_df['majority vote class']==_utils.tsg_label].index.to_series()
        novel_tsg = result_df[(result_df['majority vote class']==_utils.tsg_label) & (result_df['training list class']!=_utils.tsg_label)].index.to_series()
        pred_driver = result_df[result_df['majority vote cancer gene']==1].index.to_series()
        pred_onco.to_csv(_utils.clf_result_dir + cfg_opts['rrf_onco'], sep='\t', index=False, header=None)
        novel_onco.to_csv(_utils.clf_result_dir + cfg_opts['rrf_novel_onco'], sep='\t', index=False, header=None)
        pred_tsg.to_csv(_utils.clf_result_dir + cfg_opts['rrf_tsg'], sep='\t', index=False, header=None)
        novel_tsg.to_csv(_utils.clf_result_dir + cfg_opts['rrf_novel_tsg'], sep='\t', index=False, header=None)
        log_str = ('Majority vote Random forest: {0} ({1} novel) oncogenes, '
                   '{2} ({3} novel) tsg'.format(len(pred_onco), len(novel_onco),
                                                len(pred_tsg), len(novel_tsg)))
        logger.info(log_str)
    else:
        pred_onco = result_df[result_df['oncogene q-value']<=.1].index.to_series()
        novel_onco = result_df[(result_df['oncogene q-value']<=.1) & (result_df['training list class']!=_utils.onco_label)].index.to_series()
        pred_tsg = result_df[result_df['tsg q-value']<=.1].index.to_series()
        novel_tsg = result_df[(result_df['tsg q-value']<=.1) & (result_df['training list class']!=_utils.tsg_label)].index.to_series()
        pred_driver = result_df[result_df['driver q-value']<=.1].index.to_series()
        pred_onco.to_csv(_utils.clf_result_dir + cfg_opts['rrf_onco'], sep='\t', index=False, header=None)
        novel_onco.to_csv(_utils.clf_result_dir + cfg_opts['rrf_novel_onco'], sep='\t', index=False, header=None)
        pred_tsg.to_csv(_utils.clf_result_dir + cfg_opts['rrf_tsg'], sep='\t', index=False, header=None)
        novel_tsg.to_csv(_utils.clf_result_dir + cfg_opts['rrf_novel_tsg'], sep='\t', index=False, header=None)
        log_str = ('Random forest significance test: {0} ({1} novel) oncogenes, '
                   '{2} ({3} novel) tsg'.format(len(pred_onco), len(novel_onco),
                                                len(pred_tsg), len(novel_tsg)))
        logger.info(log_str)

    # only plot if matplotlib
    try:
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
        logger.info('Finished running sub-sampled Random Forest')

        # dummy classifier, predict most frequent
        logger.debug('Running Dummy Classifier. . .')
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
        logger.debug('Finished dummy classifier.')

        # plot oncogene roc figure
        rrandom_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_onco_mean_roc_auc
        dummy_str = 'dummy (AUC = %0.3f)' % dclf_onco_mean_roc_auc
        rrclf_onco_mean_tpr = np.mean(rrclf_onco_tpr, axis=0)
        dclf_onco_mean_tpr = np.mean(dclf_onco_tpr, axis=0)
        df = pd.DataFrame({
                        rrandom_forest_str: rrclf_onco_mean_tpr,
                        dummy_str: dclf_onco_mean_tpr},
                        index=rrclf_onco_fpr)
        line_style = {dummy_str: '--',
                    rrandom_forest_str: '-',
                    }
        save_path = _utils.clf_plot_dir + cfg_opts['roc_plot_oncogene']
        plot_data.receiver_operator_curve(df, save_path, line_style)

        # plot tsg roc figure
        r_random_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_tsg_mean_roc_auc
        dummy_str = 'dummy (AUC = %0.3f)' % dclf_tsg_mean_roc_auc
        rrclf_tsg_mean_tpr = np.mean(rrclf_tsg_tpr, axis=0)
        dclf_tsg_mean_tpr = np.mean(dclf_tsg_tpr, axis=0)
        df = pd.DataFrame({
                        r_random_forest_str: rrclf_tsg_mean_tpr,
                        dummy_str: dclf_tsg_mean_tpr},
                        index=rrclf_tsg_fpr)
        line_style = {dummy_str: '--',
                    r_random_forest_str: '-',
                    }
        save_path = _utils.clf_plot_dir + cfg_opts['roc_plot_tsg']
        plot_data.receiver_operator_curve(df, save_path, line_style)

        # plot driver roc figure
        r_random_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_driver_mean_roc_auc
        dummy_str = 'dummy (AUC = %0.3f)' % dclf_driver_mean_roc_auc
        rrclf_driver_mean_tpr = np.mean(rrclf_driver_tpr, axis=0)
        dclf_driver_mean_tpr = np.mean(dclf_driver_tpr, axis=0)
        df = pd.DataFrame({
                        r_random_forest_str: rrclf_driver_mean_tpr,
                        dummy_str: dclf_driver_mean_tpr},
                        index=rrclf_driver_fpr)
        line_style = {dummy_str: '--',
                    r_random_forest_str: '-',
                    }
        save_path = _utils.clf_plot_dir + cfg_opts['roc_plot_driver']
        plot_data.receiver_operator_curve(df, save_path, line_style)

        # plot oncogene pr figure
        rrandom_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_onco_mean_pr_auc
        dummy_str = 'dummy (AUC = %0.3f)' % dclf_onco_mean_pr_auc
        rrclf_onco_mean_precision = np.mean(rrclf_onco_precision, axis=0)
        dclf_onco_mean_precision = np.mean(dclf_onco_precision, axis=0)
        df = pd.DataFrame({
                        rrandom_forest_str: rrclf_onco_mean_precision,
                        },
                        index=rrclf_onco_recall)
        line_style = {dummy_str: '--',
                    rrandom_forest_str: '-',
                    }
        save_path = _utils.clf_plot_dir + cfg_opts['pr_plot_oncogene']
        plot_data.precision_recall_curve(df, save_path, line_style,
                                        #sem_df,
                                        title='Oncogene Precision-Recall Curve')

        # plot tsg pr figure
        r_random_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_tsg_mean_pr_auc
        dummy_str = 'dummy (AUC = %0.3f)' % dclf_tsg_mean_pr_auc
        rrclf_tsg_mean_precision = np.mean(rrclf_tsg_precision, axis=0)
        dclf_tsg_mean_precision = np.mean(dclf_tsg_precision, axis=0)
        df = pd.DataFrame({
                        r_random_forest_str: rrclf_tsg_mean_precision,
                        },
                        index=rrclf_tsg_recall)
        line_style = {dummy_str: '--',
                    r_random_forest_str: '-',
                    }
        save_path = _utils.clf_plot_dir + cfg_opts['pr_plot_tsg']
        plot_data.precision_recall_curve(df, save_path, line_style,
                                        title='TSG Precision-Recall Curve')

        # plot driver gene pr figure
        r_random_forest_str = '20/20+ Classifier (AUC = %0.3f)' % rrclf_driver_mean_pr_auc
        rrclf_driver_mean_precision = np.mean(rrclf_driver_precision, axis=0)
        df = pd.DataFrame({
                        r_random_forest_str: rrclf_driver_mean_precision,
                        },
                        index=rrclf_driver_recall)
        line_style = {dummy_str: '--',
                    r_random_forest_str: '-',
                    }
        save_path = _utils.clf_plot_dir + cfg_opts['pr_plot_driver']
        plot_data.precision_recall_curve(df, save_path, line_style,
                                        title='Driver Precision-Recall Curve')

        # save performance metrics of ROC and PR AUC
        save_path = _utils.clf_result_dir + cfg_opts['performance']
        logger.info('Saving performance metrics ({0}) . . .'.format(save_path))
        metrics = [['TSG', rrclf_tsg_mean_roc_auc, rrclf_tsg_mean_pr_auc],
                ['OG', rrclf_onco_mean_roc_auc, rrclf_onco_mean_pr_auc],
                ['Driver', rrclf_driver_mean_roc_auc, rrclf_driver_mean_pr_auc]]
        perf_df = pd.DataFrame(metrics, columns=['Type', 'ROC AUC', 'PR AUC'])
        perf_df.to_csv(save_path, sep='\t', index=False)
    except:
        pass


if __name__ == "__main__":
    main()
