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
