import utils.python.plot as myplt
import utils.python.util as _utils
import matplotlib.pyplot as plt
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def onco_mutations_parameter(df,
                             save_path,
                             title='',
                             xlabel='',
                             ylabel=''):
    logger.info('Plotting oncogene mutations while varying recurrent '
                'parameter')
    # df.index.name = 'Oncogene Score'
    myplt.line(df, save_path,
               title=title,
               ylabel=ylabel,
               xlabel=xlabel)
    logger.info('Finished plotting.')


def tsg_mutations_parameter(df,
                            save_path,
                            title='',
                            xlabel='',
                            ylabel=''):
    logger.info('Plotting tsg mutations while varying tsg score')
    myplt.line(df, save_path,
               title=title,
               ylabel=ylabel,
               xlabel=xlabel)
    logger.info('Finished plotting.')


def feature_importance_barplot(mean_df,
                               std_df,
                               save_path):
    logger.info('Plotting feature importances from random forest . . .')

    # rename columns so latex doesn't complain about '_'
    rename_dict = {col: col.replace('_', '  ') for col in mean_df.index}
    mean_df.rename(rename_dict, inplace=True)
    std_df.rename(rename_dict, inplace=True)

    mean_df.sort(ascending=True)  # sort with most important features first
    std_df = std_df.ix[mean_df.index]  # match ordering in bar plot
    myplt.barplot(mean_df,
                  save_path,
                  kind='barh',
                  xerr=std_df,
                  ecolor='#00008B',  # dark blue
                  title='Feature Importance in Random Forest',
                  xlabel='Feature Importance')
    logger.info('Finished plotting feature importance bar plot.')


def precision_recall_curve(df,
                           save_path,
                           style,
                           sem=None,
                           title='Precision-Recall Curve',
                           xlabel='Recall',
                           ylabel='Precision'):
    logger.info('Plotting precision-recall curve (%s) ...' % save_path)
    if not sem:
        # normal line plot
        myplt.line(df, save_path,
                   style=style,
                   title=title,
                   xlabel=xlabel,
                   ylabel=ylabel)
    else:
        # include standard error of the mean (sem) in the plot
        myplt.line_fill_between(df, sem, save_path,
                                style=style,
                                title=title,
                                xlabel=xlabel,
                                ylabel=ylabel)
    logger.info('Finished plotting PR curve.')


def receiver_operator_curve(df,
                            save_path,
                            style,
                            sem=None,
                            title='ROC Curve',
                            xlabel='False-Positive Rate',
                            ylabel='True-Positive Rate'):
    logger.info('Plotting receiver operator curve (%s) ...' % save_path)
    if not sem:
        # normal line plot
        myplt.line(df, save_path,
                   style=style,
                   title='ROC Curve',
                   xlabel='False-Positive Rate',
                   ylabel='True-Positive Rate')
    else:
        # include standard error of the mean (sem) in plot
        myplt.line_fill_between(df, sem, save_path,
                                title=title,
                                xlabel=xlabel,
                                ylabel=ylabel)
    logger.info('Finished plotting ROC curve.')


def vogelstein_score_scatter(df, min_count, save_path):
    df['total'] = df.sum(axis=1)
    df = df[df.total > min_count]  # filter low counts
    df['recurrent count'] = df['recurrent missense'] + df['recurrent indel']
    df['deleterious count'] = df['frame shift'] + df['nonsense'] + df['lost stop'] + df['no protein']
    df = df[(df['deleterious count']>=7) | (df['recurrent count']>=10)]
    df['oncogene score'] = df['recurrent count'].div(df['total'].astype(float))
    df['tsg score'] = df['deleterious count'].div(df['total'].astype(float))
    myplt.scatter(df['oncogene score'],
                  df['tsg score'],
                  save_path,
                  size=30,
                  title='Oncogene score vs TSG score',
                  xlabel='Oncogene score',
                  ylabel='TSG score')


def prob_kde(df, col_name, save_path,
             title,
             xlabel='Probability'):
    df['olfactory flag'] = [1 if gene in _utils.olfactory_set
                            else 0 for gene in df.index]
    try:
        df[df['olfactory flag']==1][col_name].plot(kind='kde', label='Olfactory Receptors')
    except:
        pass
    df[df['true class']==2][col_name].plot(kind='kde', label='TSG')
    df[df['true class']==1][col_name].plot(kind='kde', label='Oncogenes')
    df[df['true class']==0][col_name].plot(kind='kde', label='Other genes')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.legend()
    plt.savefig(save_path)
    plt.close()


def prob_scatter(df, plot_path, title):
    """
    Scatter plot of oncogene versus tsg probabilities.

    **Parameters**

    df : pd.DataFrame
        results from random forest prediction
    plot_path : str
        path to save scatter plot of tsg/onco probabilities
    title : str
        title of scatter plot
    """
    # scatter plot of oncogene/tsg probabilities
    myplt.scatter(df['oncogene probability'],
                  df['tsg probability'],
                  plot_path,
                  xlabel='Oncogene Probability',
                  ylabel='TSG Probability',
                  title=title,
                  colors='#348ABD')


def sample_boxplot(pred_onco,
                   pred_tsg,
                   pred_driver,
                   save_path,
                   xlabel='',
                   ylabel='',
                   title=''):
    cfg = _utils.get_output_config('sample')
    df = pd.read_csv(_utils.result_dir + cfg['max_gene_pct_sample_out'],
                     sep='\t', index_col=0)
    df['Predicted Type'] = [("oncogene" if g in pred_onco else "other") for g in df.index]
    df.ix[pred_tsg, 'Predicted Type'] = 'TSG'
    # df['tsg'] = [(1 if g in pred_tsg else 0) for g in df.index]
    # df['driver'] = [(1 if g in pred_driver else 0) for g in df.index]
    # df['other'] = 1 - df[['oncogene', 'tsg', 'driver']].max(axis=1)

    if not xlabel:
        xlabel = 'Predicted Type'
    if not ylabel:
        ylabel = 'Maximum Pct of Samples for a Tumor Type'
    if not title:
        title = 'Percentage of Samples with Non-Silent Mutation'

    myplt.boxplot(df, by='Predicted Type',
                  save_path=save_path,
                  xlabel=xlabel,
                  ylabel=ylabel,
                  title=title)
