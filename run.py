#!/usr/bin/env python
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


def start_logging(log_file='', log_level='INFO'):
    """Start logging information into the log directory.

    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """
    if not log_file:
        # create log directory if it doesn't exist
        log_dir = os.path.abspath('log') + '/'
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # specify auto-stamped log file
        log_file = 'log/log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    # logger options
    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO
    myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'

    # create logger
    if not log_file == 'stdout':
        # normal logging to a regular file
        logging.basicConfig(level=lvl,
                            format=myformat,
                            filename=log_file,
                            filemode='w')
    else:
        # logging to stdout
        root = logging.getLogger()
        root.setLevel(lvl)
        stdout_stream = logging.StreamHandler(sys.stdout)
        stdout_stream.setLevel(lvl)
        formatter = logging.Formatter(myformat)
        stdout_stream.setFormatter(formatter)
        root.addHandler(stdout_stream)
        root.propagate = True


def handle_uncaught_exceptions(t, ex, tb):
    """Handle any uncaught exceptions."""
    traceback_contents = ''.join(traceback.format_list(traceback.extract_tb(tb)))
    print('*'*40)
    print('AN ERROR HAS OCCURRED: check the log file')
    print('*'*40)
    logging.error('Type: ' + str(t))
    logging.error('Exception: ' + str(ex))
    logging.error('Traceback:\n ' + traceback_contents)
    sys.exit(EXCEPTION_EXIT_STATUS)


def _data_analysis():
    """Wrapper function to call scripts in the data_analysis folder."""
    if args.database == 'cosmic_nuc':
        # change output dir for COSMIC_nuc
        # _utils.plot_dir = 'data_analysis/plots/cosmic_nuc/'
        # _utils.result_dir = 'data_analysis/results/cosmic_nuc/'
        pass
    elif args.database == 'genes':
        # change output dir for data/genes.db
        # _utils.plot_dir = 'data_analysis/plots/genes/'
        #_utils.result_dir = 'data_analysis/results/genes/'
        pass

    src.data_analysis.python.stats.main(args.recurrent,
                                        args.recurrent_max,
                                        args.database,
                                        args.classify_only)  # run code


def _classify():
    """Wrapper function to call scripts in the classify folder."""
    opts = vars(args)  # create a dictionary for CLI options
    src.classify.python.classifier.main(opts)  # run code


def _train():
    """Wrapper function to call script in the train folder."""
    opts = vars(args)  # create a dictionary for CLI options
    src.train.python.train.main(opts)  # run code


def _simulation():
    """Wrapper function to call scripts in the classify folder."""
    opts = vars(args)  # create a dictionary for CLI options
    if not opts['random_half_split']:
        src.simulation.python.simulate_performance.main(opts)
    else:
        src.simulation.python.simulate_consistency.main(opts)


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
    parser = argparse.ArgumentParser(description='Run scripts')
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
    subparser = parser.add_subparsers(help='sub-command help')

    # data analysis sub-command
    parser_data_analysis = subparser.add_parser('data_analysis',
                                                help='Run scripts in the data'
                                                ' analysis folder')
    parser_data_analysis.set_defaults(func=_data_analysis)
    parser_data_analysis.add_argument('--classify-only',
                                      action='store_true',
                                      default=False,
                                      help='Only update information that is '
                                      'pertinent to the classify commands')
    parser_data_analysis.add_argument('-r', '--recurrent',
                                      type=int,
                                      action='store',
                                      default=2,
                                      help='Minimum number of mutations at a '
                                      'recurrent position. (default: 2)')
    parser_data_analysis.add_argument('-rm', '--recurrent-max',
                                      type=int,
                                      action='store',
                                      default=float('inf'),
                                      help='Maximum number of mutations at a '
                                      'recurrent position. (Defualt: infinity)')
    parser_data_analysis_grouper = parser_data_analysis.add_mutually_exclusive_group()
    parser_data_analysis_grouper.add_argument('-c', '--cosmic_nuc',
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
                                              ' on data/genes.db. (Default)')

    # classify sub-command
    parser_classify = subparser.add_parser('classify',
                                           help='Runs classification either with '
                                           'a provided trained classifier (using '
                                           'train) or evaluates classifier performance '
                                           'using k-fold cross-validation.')
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
                                 help='Path to file containing features in tab'
                                 'separated format. Defaults to path specified '
                                 'in config.')
    parser_classify.add_argument('-e', '--empirical-p-values',
                                 type=str,
                                 action='store',
                                 default=None,
                                 help='Path to file outputing empirical p-values '
                                 'if provided input represents simulated data. '
                                 '(Default: None)')
    parser_classify.add_argument('-s', '--simulated',
                                 action='store_true',
                                 default=False,
                                 help='Flag indicating if input features were simulated. '
                                 'Simulated data is used to construct an empirical null distribution. '
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
                                 default=3.,
                                 help='Ratio of sample size for R\'s random forest for '
                                 '"other" genes. (default: 3.0)')
    parser_classify.add_argument('-n', '--ntrees',
                                 type=int,
                                 action='store',
                                 default=200,
                                 help='Number of decision trees for random forests. '
                                 '(default: 200)')
    parser_classify.set_defaults(func=_classify)

    # train sub-command
    parser_train = subparser.add_parser('train',
                                        help='Train random forest classifier')
    parser_train.add_argument('-f', '--features',
                              type=str,
                              action='store',
                              default=None,
                              help='Path to file containing features in tab'
                                 'separated format. Defaults to path specified '
                                 'in config.')
    parser_train.add_argument('-m', '--min-count',
                              type=int,
                              action='store',
                              default=0,
                              help='Minimum number of mutations in a gene '
                              'for the gene to be considered in classification.'
                              ' (default: 0)')
    parser_train.add_argument('-d', '--driver-rate',
                              type=float,
                              action='store',
                              default=.7,
                              help='Sample rate for R\'s random forest for '
                              'oncogenes and TSGs. (default: .7)')
    parser_train.add_argument('-o', '--other-ratio',
                              type=float,
                              action='store',
                              default=3.,
                              help='Ratio of sample size for R\'s random forest for '
                              '"other" genes. (default: 3.0)')
    parser_train.add_argument('-n', '--ntrees',
                              type=int,
                              action='store',
                              default=200,
                              help='Number of decision trees for random forests. '
                              '(default: 200)')
    parser_train.add_argument('-r', '--output',
                              type=str, required=True,
                              help="Store the .Rdata file containing the trained"
                              " random forest classifier")
    parser_train.set_defaults(func=_train)

    # simulation sub-command
    parser_simulation = subparser.add_parser('simulation',
                                             help='Run simulation scripts'
                                             ' in the simulation folder')
    group = parser_simulation.add_mutually_exclusive_group(required=True)
    group.add_argument('-b', '--bootstrap',
                       action='store_true',
                       default=False,
                       help='Perform bootstrap (sample with replacement) on mutations')
    group.add_argument('-rs', '--random-samples',
                       action='store_true',
                       default=False,
                       help='Perform sample with out replacement on samples/tumors')
    group.add_argument('-rt', '--random-tumor-type',
                       action='store_true',
                       default=False,
                       help='Perform sample with out replacement based on tumor types')
    group.add_argument('-rhs', '--random-half-split',
                       action='store_true',
                       default=False,
                       help='Perform a stratified half-split of tumor samples by tumor type')
    help_str = ('If -rhs specified, use with replacement sampling if greater than zero.'
                'Else if zero, then use sampling without replacement.')
    parser_simulation.add_argument('-with-replacement', '--with-replacement',
                                   type=float,  default=0.0,
                                   help=help_str)
    parser_simulation.add_argument('-p', '--processes',
                                   type=int,
                                   action='store',
                                   default=1,
                                   help='Number of processes to use '
                                   ' for simulation (more==faster)')
    parser_simulation.add_argument('-l', '--lower-sample-rate',
                                   type=float,
                                   action='store',
                                   default=0.1,
                                   help='Lower end of sample rate interval for simulations. (Default: .1)')
    parser_simulation.add_argument('-u', '--upper-sample-rate',
                                   type=float,
                                   action='store',
                                   default=3.1,
                                   help='Upper end of sample rate interval for simulations. (Default: 3.1)')
    parser_simulation.add_argument('-num', '--num-sample-rate',
                                   type=int,
                                   action='store',
                                   default=7,
                                   help='Number of sampling rates to simulate between '
                                   'LOWER_SAMPLE_RATE and UPPER_SAMPLE_RATE. (Default: 7)')
    help_str = ('Step size for progessing further down list of top genes. Only '
                'used when random half split flag is specified')
    parser_simulation.add_argument('-step', '--step-size',
                                   type=int,
                                   default=50,
                                   help=help_str)
    help_str = ('Weight factor for ranked biased consistency measure. Only used '
                'when random half split flag is specified.')
    parser_simulation.add_argument('-weight', '--weight',
                                   type=float,
                                   default=.1,
                                   help=help_str)
    help_str = ('Maximum depth of genes from top of list to consider for consistency. '
                'Only used when random half split flag is specified.')
    parser_simulation.add_argument('-depth', '--depth',
                                   type=int,
                                   default=200,
                                   help=help_str)
    parser_simulation.add_argument('-m', '--min-count',
                                   type=int,
                                   action='store',
                                   default=0,
                                   help='Minimum number of mutations in a gene '
                                   'for the gene to be considered in classification.'
                                   ' (default: 0)')
    parser_simulation.add_argument('-d', '--driver-rate',
                                   type=float,
                                   action='store',
                                   default=.7,
                                   help='Sample rate for R\'s random forest for '
                                   'oncogenes and TSGs. (default: .7)')
    parser_simulation.add_argument('-o', '--other-ratio',
                                   type=float,
                                   action='store',
                                   default=3.,
                                   help='Ratio of sample size for R\'s random forest for '
                                   '"other" genes. (default: 3.0)')
    parser_simulation.add_argument('-n', '--ntrees',
                                   type=int,
                                   action='store',
                                   default=500,
                                   help='Number of decision trees for random forests. '
                                   '(default: 500)')
    parser_simulation.add_argument('-s', '--samples',
                                   type=int,
                                   action='store',
                                   default=10,
                                   help='Number of samples for each simulation.')
    parser_simulation.add_argument('--betweeness',
                                   action='store_true',
                                   default=False,
                                   help='Gene betweeness from Biogrid.')
    parser_simulation.add_argument('--degree',
                                   action='store_true',
                                   default=False,
                                   help='Gene edge degree.')
    parser_simulation.add_argument('--gene-length',
                                   action='store_true',
                                   default=False,
                                   help='Add gene length to features for '
                                   'simulation command')
    parser_simulation.add_argument('--mutation-rate',
                                   action='store_true',
                                   default=False,
                                   help='Add noncoding mutation rate to'
                                   ' features for simulation command')
    parser_simulation.add_argument('--replication-time',
                                   action='store_true',
                                   default=False,
                                   help='Add replication time to'
                                   ' features for simulation command')
    parser_simulation.add_argument('--expression',
                                   action='store_true',
                                   default=False,
                                   help='Add gene expression to'
                                   ' features for simulation command')
    parser_simulation.set_defaults(func=_simulation)

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
                                           help=help_string)
    parser_features.set_defaults(func=_features)
    parser_features.add_argument('-m', '--min-count',
                                 type=int,
                                 action='store',
                                 default=0,
                                 help='Minimum number of mutations in a gene '
                                 'for the gene to be in the saved feature file.'
                                 ' (default: 0)')
    parser_features.add_argument('--betweeness',
                                 action='store_true',
                                 default=False,
                                 help='Gene betweeness from Biogrid.')
    parser_features.add_argument('--degree',
                                 action='store_true',
                                 default=False,
                                 help='Gene edge degree.')
    parser_features.add_argument('--gene-length',
                                 action='store_true',
                                 default=False,
                                 help='Add gene length to features for '
                                 'classify command')
    parser_features.add_argument('--mutation-rate',
                                 action='store_true',
                                 default=False,
                                 help='Add noncoding mutation rate to'
                                 ' features for classify command')
    parser_features.add_argument('--replication-time',
                                 action='store_true',
                                 default=False,
                                 help='Add replication time to'
                                 ' features for classify command')
    parser_features.add_argument('--expression',
                                 action='store_true',
                                 default=False,
                                 help='Add gene expression to'
                                 ' features for classify command')

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
    log_level = args.log_level
    start_logging(log_file=log_file,
                  log_level=log_level)  # start logging

    # import all the modules for 20/20+
    import src.data_analysis.python.stats
    import src.classify.python.classifier
    import src.simulation.python.simulate_performance
    import src.simulation.python.simulate_consistency
    import src.features.python.features
    import src.savedb.python.gene_tsv
    import src.savedb.python.gene_features
    import src.savedb.python.gene_maf
    import src.savedb.python.merge_mutations
    import src.train.python.train
    import src.utils.python.util as _utils

    # make output directory if specified by user
    save_dir = args.out_dir
    if save_dir is not None:
        _opts = _utils.get_input_config('result')
        _utils.plot_dir = os.path.join(save_dir, _opts['plot_dir'])
        _utils.result_dir = os.path.join(save_dir, _opts['result_dir'])
        _utils.clf_plot_dir = os.path.join(save_dir, _opts['clf_plot_dir'])
        _utils.clf_result_dir = os.path.join(save_dir, _opts['clf_result_dir'])
        _utils.feature_plot_dir = os.path.join(save_dir, _opts['feature_plot_dir'])
        _utils.sim_plot_dir = os.path.join(save_dir, _opts['sim_plot_dir'])
        _utils.sim_result_dir = os.path.join(save_dir, _opts['sim_result_dir'])
        if not os.path.exists(_utils.plot_dir): os.makedirs(_utils.plot_dir)
        if not os.path.exists(_utils.result_dir): os.makedirs(_utils.result_dir)
        if not os.path.exists(_utils.clf_plot_dir): os.makedirs(_utils.clf_plot_dir)
        if not os.path.exists(_utils.clf_result_dir): os.makedirs(_utils.clf_result_dir)
        if not os.path.exists(_utils.feature_plot_dir): os.makedirs(_utils.feature_plot_dir)
        if not os.path.exists(_utils.sim_plot_dir): os.makedirs(_utils.sim_plot_dir)
        if not os.path.exists(_utils.sim_result_dir): os.makedirs(_utils.sim_result_dir)

    args.func()  # run function corresponding to user's command
    logging.info('FINISHED SUCCESSFULLY!')
