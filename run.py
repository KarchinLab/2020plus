import logging
import datetime
import os
import sys
import traceback
import argparse
import data_analysis.python.stats

# define exit status
EXCEPTION_EXIT_STATUS = 1
BAD_ARG_EXIT_STATUS = 2


def start_logging(log_file=''):
    """Start logging information into the log directory.

    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """
    if not log_file:
        log_file = 'log/log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    logging.basicConfig(level=logging.INFO,
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
    if args.stats:
        data_analysis.python.stats.main()


if __name__ == '__main__':
    # initializations
    sys.excepthook = handle_uncaught_exceptions  # handle exceptions

    parser = argparse.ArgumentParser(description='Run scripts')
    parser.add_argument('-l', '--log',
                        action='store_true',
                        help='write a log file')
    subparser = parser.add_subparsers(help='sub-command help')
    parser_data_analysis = subparser.add_parser('data_analysis',
                                                help='Run scripts in data'
                                                'analysis folder')
    parser_data_analysis.set_defaults(func=_data_analysis)
    parser_data_analysis.add_argument('-s', '--stats',
                                      action='store_true',
                                      help='Generate data analysis stats')

    args = parser.parse_args()
    log_file = '' if args.log else os.devnull
    start_logging(log_file=log_file)  # start logging
    args.func()
    logging.info('FINISHED SUCCESSFULLY!')


    # run program
    # logger = logging.getLogger(__name__)
