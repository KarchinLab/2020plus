import logging
import datetime
import os
import sys
import traceback
import argparse
import data_analysis.python.stats
import classify.python.classifier
import utils.python.gene_tsv
import utils.python.util as _utils

# define exit status
EXCEPTION_EXIT_STATUS = 1
BAD_ARG_EXIT_STATUS = 2


def start_logging(log_file='', log_level='INFO'):
    """Start logging information into the log directory.

    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """
    if not log_file:
        log_file = 'log/log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO
    logging.basicConfig(level=lvl,
                        format='%(asctime)s - %(name)s - '
                               '%(levelname)s - %(message)s',
                        filename=log_file,
                        filemode='w')


def handle_uncaught_exceptions(t, ex, tb):
    """Handle any uncaught exceptions."""
    traceback_contents = ''.join(traceback.format_list(traceback.extract_tb(tb)))
    print(t)
    print(traceback_contents)
    print(ex)
    logging.error('Type: ' + str(t))
    logging.error('Exception: ' + str(ex))
    logging.error('Traceback:\n ' + traceback_contents)
    sys.exit(EXCEPTION_EXIT_STATUS)


def _data_analysis():
    """Wrapper function to call scripts in the data_analysis folder."""
    if args.database == 'cosmic_nuc':
        # change output dir for COSMIC_nuc
        _utils.plot_dir = 'data_analysis/plots/cosmic_nuc/'
        _utils.result_dir = 'data_analysis/results/cosmic_nuc/'
    elif args.database == 'genes':
        # change output dir for data/genes.db
        _utils.plot_dir = 'data_analysis/plots/genes/'
        _utils.result_dir = 'data_analysis/results/genes/'

    data_analysis.python.stats.main(args.database)  # run code

def _classify():
    """Wrapper function to call scripts in the classify folder."""
    classify.python.classifier.main()


def _savedb():
    """Wrapper function to call gene_tsv's main function"""
    utils.python.gene_tsv.main()


if __name__ == '__main__':
    # initializations
    sys.excepthook = handle_uncaught_exceptions  # handle exceptions

    parser = argparse.ArgumentParser(description='Run scripts')
    parser.add_argument('-l', '--log',
                        type=str,
                        action='store',
                        default='',
                        help='Write a log file (--log=DEBUG for debug mode, '
                        '--log=INFO for info mode)')
    subparser = parser.add_subparsers(help='sub-command help')
    parser_data_analysis = subparser.add_parser('data_analysis',
                                                help='Run scripts in the data'
                                                ' analysis folder')
    parser_data_analysis.set_defaults(func=_data_analysis)
    parser_data_analysis_grouper = parser_data_analysis.add_mutually_exclusive_group()
    parser_data_analysis_grouper.add_argument('-cn', '--cosmic_nuc',
                                              dest='database',
                                              action='store_const',
                                              const='cosmic_nuc',
                                              help='Generate data analysis stats'
                                              ' on COSMIC_nuc database')
    parser_data_analysis_grouper.add_argument('-g', '--genes',
                                              dest='database',
                                              action='store_const',
                                              const='genes',
                                              help='Generate data analysis stats'
                                              ' on data/genes.db')
    parser_classify = subparser.add_parser('classify',
                                           help='Run classification scripts'
                                           ' in the classify folder')
    parser_classify.set_defaults(func=_classify)
    help_string = ('Concatenate tab delim gene files found in /databases/COSMIC '
                   'and then save them to a sqlite database for further use.')
    parser_savedb = subparser.add_parser('savedb',
                                         help=help_string)
    parser_savedb.set_defaults(func=_savedb)

    parser.set_defaults(database='genes')  # by default work on sqlite db

    args = parser.parse_args()
    log_file = '' if args.log else os.devnull
    log_level = args.log
    start_logging(log_file=log_file,
                  log_level=log_level)  # start logging
    args.func()
    logging.info('FINISHED SUCCESSFULLY!')
