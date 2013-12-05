import utils.python.plot as myplt
import logging

logger = logging.getLogger(__name__)

def onco_mutations_parameter(df,
                             save_path,
                             title='',
                             xlabel='',
                             ylabel=''):
    logger.info('Plotting oncogene mutations while varying recurrent '
                'parameter')
    myplt.line(df, save_path,
               title=title,
               ylabel=ylabel,
               xlabel=xlabel)
    logger.info('Finished plotting.')


def feature_importance_barplot(mean_df,
                               std_df,
                               save_path):
    logger.info('Plotting feature importances from random forest . . .')
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


def precision_recall_curve(df, save_path, style,
                           title='Precision-Recall Curve',
                           xlabel='Recall',
                           ylabel='Precision'):
    logger.info('Plotting precision-recall curve (%s) ...' % save_path)
    myplt.line(df, save_path,
               style=style,
               title=title,
               xlabel=xlabel,
               ylabel=ylabel)
    logger.info('Finished plotting PR curve.')


def receiver_operator_curve(df, save_path, style,
                            title='ROC Curve',
                            xlabel='False-Positive Rate',
                            ylabel='True-Positive Rate'):
    logger.info('Plotting receiver operator curve (%s) ...' % save_path)
    myplt.line(df, save_path,
               style=style,
               title='ROC Curve',
               xlabel='False-Positive Rate',
               ylabel='True-Positive Rate')
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

