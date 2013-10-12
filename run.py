import logging
import datetime
import os
import sys
import traceback

# define exit status
EXCEPTION_EXIT_STATUS = 1
BAD_ARG_EXIT_STATUS = 2

def start_logging():
    """Start logging information into the log directory."""
    log_file = 'log/log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
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

if __name__=='__main__':
    # initializations
    start_logging()  # start logging
    sys.excepthook = handle_uncaught_exceptions  # handle exceptions

    raise ValueError('An error occurred')

    # run program
    logger = logging.getLogger(__name__)
    logger.info('what is up?')
    logger.debug('there is a bug')
    try:
        raise ValueError
    except:
        logger.error('an error occured', exc_info=True)
