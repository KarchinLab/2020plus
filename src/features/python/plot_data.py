import src.utils.python.util as _utils
import src.utils.python.plot as myplt
import logging

logger = logging.getLogger(__name__)

def or_gene_hist(feature_df, save_path,
                 title='Mutations for Olfactory Receptor Genes',
                 xlabel='Number of Mutations',
                 ylabel='Number of Genes'):
    logger.info('Plotting Olfactory Receptor histogram (%s) . . .' % save_path)
    feature_df = feature_df.set_index('gene')  # make gene name the index
    or_df = feature_df.ix[_utils.olfactory_set]  # only olfactory genes
    myplt.histogram(or_df['total'], save_path,
                    title=title,
                    xlabel=xlabel,
                    ylabel=ylabel)
    logger.info('Finsihed plotting.')


def correlation_plot(feature_df,
                     column_x,
                     column_y,
                     save_path,
                     title='',
                     xlabel='',
                     ylabel=''):
    logger.info('Plotting correlation between {0} and {1}.'.format(column_x, column_y))
    myplt.correlation_plot(feature_df[column_x],
                           feature_df[column_y],
                           save_path,
                           title=title,
                           xlabel=xlabel,
                           ylabel=ylabel)
    #myplt.scatter(feature_df[column_x],
                  #feature_df[column_y],
                  #save_path,
                  #title=title,
                  #xlabel=xlabel,
                  #ylabel=ylabel)
    logger.info('Finished plotting correlation between '
                '{0} and {1}.'.format(column_x, column_y))
