import utils.python.util as _utils
import utils.python.plot as myplt
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
