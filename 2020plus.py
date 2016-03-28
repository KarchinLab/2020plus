#!/usr/bin/env python
# import print function for printing to stderr
from __future__ import print_function
# force project root directory to be in path. Otherwise
# package imports will fail if run.py is ran from another
# directory.
import os
import sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

# regular imports
import logging
import datetime
import traceback
import argparse

# define exit status
EXCEPTION_EXIT_STATUS = 1
BAD_ARG_EXIT_STATUS = 2


def handle_uncaught_exceptions(t, ex, tb):
    """Handle any uncaught exceptions."""
    traceback_contents = ''.join(traceback.format_list(traceback.extract_tb(tb)))
    print('*'*40, file=sys.stderr)
    print('AN ERROR HAS OCCURRED: check the log file', file=sys.stderr)
    print('*'*40, file=sys.stderr)
    logging.error('Type: ' + str(t))
    logging.error('Exception: ' + str(ex))
    logging.error('Traceback:\n ' + traceback_contents)
    sys.exit(EXCEPTION_EXIT_STATUS)


def _classify():
    """Wrapper function to call scripts in the classify folder."""
    opts = vars(args)  # create a dictionary for CLI options
    src.classify.python.classifier.main(opts)  # run code


def _train():
    """Wrapper function to call script in the train folder."""
    opts = vars(args)  # create a dictionary for CLI options
    src.train.python.train.main(opts)  # run code


def _savedb():
    """Wrapper function to call gene_tsv/gene_features main function.

    Saves information mostly from COSMIC into a database. Additional
    information from the MutSigCV paper is also stored in the gene_features
    table.
    """
    src.savedb.python.gene_tsv.main(args.hypermutator,
                                    # args.cell_line,
                                    args.input,
                                    args.output,
                                    args.no_cosmic,
                                    vars(args))  # populate the nucleotide table
    src.savedb.python.gene_features.main(args.output)  # populate the gene_features table
    src.savedb.python.gene_maf.main(args.maf,
                                    args.output,
                                    args.hypermutator)
    src.savedb.python.merge_mutations.main(args.output)


def _features():
    """Wrapper function to call the features main function."""
    opts = vars(args)  # make CLI options a dictionary
    src.features.python.features.main(opts)


if __name__ == '__main__':
    # initializations
    sys.excepthook = handle_uncaught_exceptions  # handle exceptions

    # parse command line arguments
    parser = argparse.ArgumentParser(description='Run 20/20+ pipeline')
    parser.add_argument('--out-dir',
                        type=str,
                        action='store',
                        default=None,
                        help='Path to output directory. Used by all positional arguments. '
                        '(Default: result/)')
    parser.add_argument('-ll', '--log-level',
                        type=str,
                        action='store',
                        default='',
                        help='Write a log file (--log-level=DEBUG for debug mode, '
                        '--log-level=INFO for info mode)')
    parser.add_argument('-l', '--log',
                        type=str,
                        action='store',
                        default='',
                        help='Path to log file. (accepts stdout)')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False,
                        help='Flag for more verbose log output')
    subparser = parser.add_subparsers(help='sub-command help')

    # savedb sub-command
    help_string = ('Concatenate tab delim gene files found in /databases/COSMIC '
                   'and then save them to a sqlite database for further use. '
                   'Gene length and information from the MutSigCV paper are also '
                   'stored in the database.')
    parser_savedb = subparser.add_parser('savedb',
                                         help=help_string)
    parser_savedb.set_defaults(func=_savedb)
    help_str = ('Save sqlite db to non-standard location which'
                ' will NOT be used for further analysis (optional)')
    parser_savedb.add_argument('-o', '--output',
                               type=str, default='',
                               help=help_str)
    help_str = ('Mutations from COSMIC. Either a directory unpacked from genes.tgz file or '
                'CosmicMutantExport.tsv (optional)')
    parser_savedb.add_argument('-i', '--input',
                               type=str, default='',
                               help=help_str)
    help_str = 'Flag indicating whether to only use genome wide screens in COSMIC'
    parser_savedb.add_argument('--only-genome-wide',
                               action='store_true',
                               default=False,
                               help=help_str)
    help_str = 'Include variants with unknown somatic status (potentially germline)'
    parser_savedb.add_argument('--use-unknown-status',
                               action='store_true',
                               default=False,
                               help=help_str)
    parser_savedb.add_argument('-hm', '--hypermutator',
                               type=int,
                               action='store',
                               default=500,
                               help='Number of mutations to define a sample '
                               'as a hypermutator. Hypermutator samples are filtered '
                               ' from further analysis. (default: 500)')
    parser_savedb.add_argument('-m', '--maf',
                               type=str, default='',
                               help='MAF file to augment mutations '
                               'found in COSMIC')
    parser_savedb.add_argument('-nc', '--no-cosmic',
                               action='store_true',
                               default=False,
                               help='Don\'t use mutations from COSMIC')

    # features sub-command
    help_string = ('Generate the features used in classification.'
                   ' This command should be ran before "classify".'
                   ' Features are saved as a text file.')
    parser_features = subparser.add_parser('features',
                                           help=help_string,
                                           description=help_string)
    parser_features.set_defaults(func=_features)
    help_str = 'mutation annotate output from probabilistic 20/20'
    parser_features.add_argument('-s', '--summary',
                        type=str, required=True,
                        help=help_str)
    help_str = 'TSG output from probabilistic 20/20 ("probabilistic2020 tsg")'
    parser_features.add_argument('-tsg-test', '--tsg-test',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Oncogene output from probabilistic 20/20 ("probabilistic2020 oncogene")'
    parser_features.add_argument('-og-test', '--og-test',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Mutsigcv covariate features (Default: use config file)'
    parser_features.add_argument('-c', '--covariates',
                        type=str, default=None,
                        help=help_str)
    help_str = 'BioGrid interaction network statistics (Default: use config file)'
    parser_features.add_argument('-b', '--biogrid',
                        type=str, default=None,
                        help=help_str)
    help_str = 'Randomly permute biogrid features (use for null distribution only)'
    parser_features.add_argument('-p', '--permute-biogrid',
                        action='store_true', default=False,
                        help=help_str)
    parser_features.add_argument('-rs', '--random-seed',
                                 type=int, action='store',
                                 default=71,
                                 help='Random seed for permute biogrid option (default: 71)')
    help_str = 'Output feature file for 20/20+'
    parser_features.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)

    # train sub-command
    parser_train = subparser.add_parser('train',
                                        help='Train random forest classifier (only used for null distribution)',
                                        description='Train random forest classifier (only used for null distribution)')
    parser_train.add_argument('-f', '--features',
                              type=str,
                              action='store', required=True,
                              help='Path to file containing features in tab '
                              'separated format. Defaults to path specified '
                              'in config.')
    parser_train.add_argument('-d', '--driver-rate',
                              type=float,
                              action='store',
                              default=.7,
                              help='Sample rate for R\'s random forest for '
                              'oncogenes and TSGs. (default: .7)')
    parser_train.add_argument('-o', '--other-ratio',
                              type=float,
                              action='store',
                              default=1.,
                              help='Ratio of sample size for R\'s random forest for '
                              '"other" genes. (default: 1.0)')
    parser_train.add_argument('-n', '--ntrees',
                              type=int,
                              action='store',
                              default=500,
                              help='Number of decision trees for random forests. '
                              '(default: 500)')
    parser_train.add_argument('-m', '--min-count',
                              type=int,
                              action='store',
                              default=0,
                              help='Minimum number of mutations in a gene '
                              'for the gene to be considered in classification.'
                              ' (default: 0)')
    parser_train.add_argument('-rs', '--random-seed',
                              type=int, action='store',
                              default=71,
                              help='Random seed (default: 71)')
    parser_train.add_argument('-r', '--output',
                              type=str, required=True,
                              help="Store the .Rdata file containing the trained"
                              " random forest classifier")
    parser_train.set_defaults(func=_train)

    # classify sub-command
    parser_classify = subparser.add_parser('classify',
                                           help='Runs classification either with '
                                           'a provided trained classifier (using '
                                           'train) or '
                                           'using k-fold cross-validation (no train command needed).')
    parser_classify.add_argument('-t', '--trained-classifier',
                                 type=str,
                                 action='store',
                                 default=None,
                                 help='If provided, use trained classifier from '
                                 'the train sub-command. Otherwise, perform '
                                 'cross-validation within the data set (default: None)')
    parser_classify.add_argument('-f', '--features',
                                 type=str,
                                 action='store',
                                 default=None,
                                 help='Path to file containing features in tab '
                                 'separated format. Defaults to path specified '
                                 'in config.')
    parser_classify.add_argument('-nd', '--null-distribution',
                                 type=str,
                                 action='store',
                                 default=None,
                                 help='Path to file outputing null distiribution for p-values '
                                 'if provided input represents simulated data. '
                                 '(Default: None)')
    parser_classify.add_argument('-s', '--simulated',
                                 action='store_true',
                                 default=False,
                                 help='Flag indicating if input features were simulated. '
                                 'Simulated data is used to construct a null distribution. '
                                 '(Default: False)')
    parser_classify.add_argument('-m', '--min-count',
                                 type=int,
                                 action='store',
                                 default=0,
                                 help='Minimum number of mutations in a gene '
                                 'for the gene to be considered in classification.'
                                 ' (default: 0)')
    parser_classify.add_argument('-d', '--driver-rate',
                                 type=float,
                                 action='store',
                                 default=.7,
                                 help='Sample rate for R\'s random forest for '
                                 'oncogenes and TSGs. (default: .7)')
    parser_classify.add_argument('-o', '--other-ratio',
                                 type=float,
                                 action='store',
                                 default=1.,
                                 help='Ratio of sample size for R\'s random forest for '
                                 '"other" genes. (default: 1.0)')
    parser_classify.add_argument('-n', '--ntrees',
                                 type=int,
                                 action='store',
                                 default=200,
                                 help='Number of decision trees for random forests. '
                                 '(default: 200)')
    parser_classify.add_argument('-rs', '--random-seed',
                                 type=int, action='store',
                                 default=71,
                                 help='Random seed (default: 71)')
    parser_classify.set_defaults(func=_classify)

    parser.set_defaults(database='genes')  # by default work on sqlite db
    args = parser.parse_args()  # parse the command line options

    # handle logging
    if args.log_level or args.log:
        if args.log:
            log_file = args.log
        else:
            log_file = ''  # auto-name the log file
    else:
        log_file = os.devnull
    import src.utils.python.util as _utils
    log_level = args.log_level
    _utils.start_logging(log_file=log_file,
                         log_level=log_level,
                         verbose=args.verbose)  # start logging

    # log user entered command
    logging.info('Command: {0}'.format(' '.join(sys.argv)))

    # import all the modules for 20/20+
    import src.classify.python.classifier
    import src.features.python.features
    import src.savedb.python.gene_tsv
    import src.savedb.python.gene_features
    import src.savedb.python.gene_maf
    import src.savedb.python.merge_mutations
    import src.train.python.train

    # make output directory if specified by user
    save_dir = args.out_dir
    _utils.make_result_dir(save_dir)

    args.func()  # run function corresponding to user's command
    logging.info('FINISHED SUCCESSFULLY!')
