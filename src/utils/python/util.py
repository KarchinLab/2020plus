import numpy as np
import pandas as pd
from amino_acid import AminoAcid
from nucleotide import Nucleotide
from cosmic_db import get_cosmic_db
import sqlite3
import pandas.io.sql as psql
import ConfigParser
import logging
import os
from functools import wraps

logger = logging.getLogger(__name__)

proj_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..'))
config_dir = os.path.join(proj_dir, 'config/')

onco_label, tsg_label, other_label = 1, 2, 0
class_to_label = {'oncogene': onco_label,
                  'tsg': tsg_label,
                  'other': other_label}

def get_input_config(section):
    """Returns the config object to input.cfg."""
    cfg = ConfigParser.ConfigParser()
    cfg.read(config_dir + 'input.cfg')
    cfg_options = dict(cfg.items(section))
    return cfg_options


# setup directory paths
_opts = get_input_config('result')
save_dir = os.path.join(proj_dir, _opts['save_dir'])
plot_dir = save_dir + _opts['plot_dir']
result_dir = save_dir + _opts['result_dir']
clf_plot_dir = save_dir + _opts['clf_plot_dir']
clf_result_dir = save_dir + _opts['clf_result_dir']
feature_plot_dir = save_dir + _opts['feature_plot_dir']
sim_plot_dir = save_dir + _opts['sim_plot_dir']
sim_result_dir = save_dir + _opts['sim_result_dir']

def start_logging(log_file='', log_level='INFO', verbose=False):
    """Start logging information into the log directory.

    If os.devnull is specified as the log_file then the log file will
    not actually be written to a file.
    """
    if not log_file:
        # create log directory if it doesn't exist
        log_dir = os.path.abspath('log') + '/'
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)

        # path to new log file
        log_file = log_dir + 'log.run.' + str(datetime.datetime.now()).replace(':', '.') + '.txt'

    # logger options
    lvl = logging.DEBUG if log_level.upper() == 'DEBUG' else logging.INFO

    # define logging format
    if verbose:
        myformat = '%(asctime)s - %(name)s - %(levelname)s \n>>>  %(message)s'
    else:
        myformat = '%(message)s'

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


def log_error_decorator(f):
    """Writes exception to log file if occured in decorated function.

    This decorator wrapper is needed for multiprocess logging since otherwise
    the python multiprocessing module will obscure the actual line of the error.
    """
    @wraps(f)
    def wrapper(*args, **kwds):
        try:
            result = f(*args, **kwds)
            return result
        except KeyboardInterrupt:
            logger.info('Ctrl-C stopped a process.')
        except Exception, e:
            logger.exception(e)
            raise
    return wrapper


def keyboard_exit_wrapper(func):
    def wrap(self, timeout=None):
        # Note: the timeout of 1 googol seconds introduces a rather subtle
        # bug for Python scripts intended to run many times the age of the universe.
        return func(self, timeout=timeout if timeout is not None else 1e100)
    return wrap


def read_aa_properties(file_path):
    """Read aa property counts from the data_analysis/results folder.

    Parameters
    ----------
    file_path : str
        path to aa_change.properties.txt

    Returns
    -------
    df : pd.DataFrame
        contains mutation counts for amino acid chemical properties
    """
    logger = logging.getLogger(name=__name__)
    logger.info('reading in %s ...' % file_path)
    df = pd.read_csv(file_path, sep='\t')  # read file
    df = df.set_index('initial_prop')  # set rows as initial property
    logger.info('finished reading file.')
    return df


def read_oncogenes():
    """Reads in the oncogenes from vogelsteins' science paper.

    Oncogenes from supplementary 2A:
    http://www.sciencemag.org/content/339/6127/1546.full

    Returns
    -------
    oncogenes : tuple
        tuple of gene names considered oncogenes
    """
    cfg_opts = get_input_config('input')
    onco_path = os.path.join(proj_dir, cfg_opts['oncogene'])
    with open(onco_path, 'r') as handle:
        oncogenes = tuple(gene.strip() for gene in handle.readlines())
    return oncogenes


def read_tsgs():
    """Reads in the tumor suppressor genes from vogelsteins' science paper.

    TSGs from supplementary 2A:
    http://www.sciencemag.org/content/339/6127/1546.full

    Returns
    -------
    tsgs : tuple
        tuple of gene names considered as tumor suppressors
    """
    cfg_opts = get_input_config('input')
    tsg_path = os.path.join(proj_dir, cfg_opts['tsg'])
    with open(tsg_path, 'r') as handle:
        tsgs = tuple(gene.strip() for gene in handle.readlines())
    return tsgs


def read_smg():
    """Reads in the significantly mutated genes from kandoth et al's
    nature paper.

    The paper was from the Pan-Cancer effort from TCGA.

    SMGs from supplementary:
    http://www.nature.com/nature/journal/v502/n7471/full/nature12634.html

    Returns
    -------
    smgs : tuple
        tuple of gene names considered as significantly mutated
    """
    cfg_opts = get_input_config('input')
    smg_path = os.path.join(proj_dir, cfg_opts['smg'])
    with open(smg_path, 'r') as handle:
        smgs = tuple(gene.strip() for gene in handle.readlines())

    # open connection
    try:
        # if DB is not created this will throw an error
        gene_db_path = os.path.join(proj_dir, get_db_config('2020plus')['db'])
        conn = sqlite3.connect(gene_db_path)

        sql = ("SELECT DISTINCT Gene"
              " FROM mutation"
              " WHERE Gene in " + str(smgs))
        df = psql.frame_query(sql, con=conn)
        conn.close()  # close connection

        # get significantly mutated genes found in database
        smgs_in_database = tuple(df['Gene'].astype(str))
        logger.debug('There are only %d/%d significantly mutated genes found in the database.'
                     % (len(smgs_in_database), len(smgs)))
    except:
        smgs_in_database = smgs
    return smgs_in_database


def read_cgc():
    """Gets the genes from the cancer gene census.

    Data from CGC is available from here:
    http://cancer.sanger.ac.uk/cancergenome/projects/census/

    Returns
    -------
    cgc_in_database : tuple
        tuple of gene names in cancer gene census also in COSMIC
    """
    cfg_opts = get_input_config('input')
    cgc_path = os.path.join(proj_dir, cfg_opts['cgc'])
    with open(cgc_path, 'r') as handle:
        cgc = tuple(gene.strip() for gene in handle.readlines())

    # open connection
    try:
        # if DB is not created this will throw an error
        gene_db_path = os.path.join(proj_dir, get_db_config('2020plus')['db'])
        conn = sqlite3.connect(gene_db_path)

        sql = ("SELECT DISTINCT Gene"
              " FROM mutation"
              " WHERE Gene in " + str(cgc))
        df = psql.frame_query(sql, con=conn)
        conn.close()  # close connection

        # get significantly mutated genes found in database
        cgc_in_database = tuple(df['Gene'].astype(str))
        logger.debug('There are only %d/%d CGC genes found in the database.'
                     % (len(cgc_in_database), len(cgc)))
    except:
        cgc_in_database = cgc
    return cgc_in_database


def read_olfactory_receptors():
    """Reads in the significant olfactory receptors from Mutsigcv.

    Returns
    -------
    olfactory : tuple
        tuple of gene names considered as olfactory receptors
    """
    cfg_opts = get_input_config('input')
    or_path = os.path.join(proj_dir, cfg_opts['olfactory_receptors'])
    with open(or_path, 'r') as handle:
        olfactory = tuple(gene.strip() for gene in handle.readlines())

    # open connection
    try:
        # it table is not found just catch exception
        gene_db_path = os.path.join(proj_dir, get_db_config('2020plus')['db'])
        conn = sqlite3.connect(gene_db_path)

        sql = ("SELECT DISTINCT Gene"
               " FROM mutation"
               " WHERE Gene in " + str(olfactory))

        df = psql.frame_query(sql, con=conn)
        conn.close()  # close connection

        # get significantly mutated genes found in database
        olfactory_in_database = tuple(df['Gene'])
        logger.debug('There are only %d/%d olfactory receptors found in the database.'
                     % (len(olfactory_in_database), len(olfactory)))
    except:
        olfactory_in_database = olfactory
    return olfactory_in_database


def classify_gene(gene):
    """Return whether the gene is an oncogene, tsg, or other.

    Parameters
    ----------
    gene : str
        Official gene name

    Returns
    -------
        Str, ['oncogene' | 'tsg' | 'other']
    """
    if gene in oncogene_set:
        return 'oncogene'
    elif gene in tsg_set:
        return 'tsg'
    else:
        return 'other'


def get_mutation_types(mut_iterable,
                       dna_series=None,
                       known_type=None,
                       kind='amino acid'):
    """Classify each protein HGVS mutation as a certain type.

    Parameters
    ----------
    mut_iterable : iterable
        iterable container with HGVS mutaiton strings. If amino acids,
        a secon dna_iterable is needed to identify splice mutations.
    dna_series : pd.Series
        optional, only required to find splice mutations when a list of
        amino acids is given for mut_iterable
    known_type : pd.Series
        contains list of mutation types

    Returns
    -------
    mut_type_series : pd.Series
        container of protein mutation types in same order as input
    """
    mut_type = []
    if kind == 'amino acid':
        if dna_series is None:
            # dna iterable required
            raise ValueError('DNA should be specified to identify splice mutations.')
        for i, hgvs_aa in enumerate(mut_iterable):
            aa = AminoAcid(hgvs=hgvs_aa)
            nuc = Nucleotide(hgvs=dna_series.iloc[i])
            if nuc.is_splicing_mutation:
                # check if mutation in splice site
                mut_type.append('Splice_Site')
            elif known_type is not None and known_type.iloc[i]=='Splice_Site':
                mut_type.append('Splice_Site')
            else:
                # if not in splice site, just add
                mut_type.append(aa.mutation_type)
    elif kind == 'nucleotide':
        for hgvs_nuc in mut_iterable:
            nuc = Nucleotide(hgvs=hgvs_nuc)
            mut_type.append(nuc.mutation_type)
    mut_type_series = pd.Series(mut_type)
    return mut_type_series


def count_mutation_types(hgvs_iterable, dna_series=None, known_type=None, kind='amino acid'):
    """Count mutation types from HGVS protein strings (missense, indels, etc.)
    and DNA strings (substitutions, indels).

    Parameters
    ----------
    hgvs_iterable : iterable
        An iterable object containing protein HGVS
    dna_iterable : iterable
        contains hgvs DNA mutations to classify splice mutations
        for amino acid. Only required if hgvs_iterable is AA mutations.
    known_type : pd.Series
        known mutation consequence type

    Returns
    -------
    unique_cts : pd.Series
        A pandas series object counting protein mutation types
    """
    mut_type_series = get_mutation_types(hgvs_iterable,
                                         dna_series=dna_series,
                                         known_type=known_type,
                                         kind=kind)  # get mutation types
    unique_cts = mut_type_series.value_counts() # count mutation types
    return unique_cts


def get_output_config(section):
    """Returns the config object to output.cfg."""
    cfg = ConfigParser.ConfigParser()
    cfg.read(config_dir + 'output.cfg')
    cfg_options = dict(cfg.items(section))
    return cfg_options


def get_db_config(section):
    """Return the config object to db.cfg."""
    cfg = ConfigParser.ConfigParser()
    cfg.read(config_dir + 'db.cfg')
    cfg_options = dict(cfg.items(section))
    return cfg_options


def read_cosmic_tsv_by_gene(gene_name):
    """Reads the stored flat file corresponding to the gene_name.

    NOTE: Assumes cosmic flat files are in cosmic_dir specified by input.cfg
    and are sorted into alphabetical directories (eg. 'A'...'Z').

    Parameters
    ----------
    gene_name : str
        gene name

    Returns
    -------
    df : pd.DataFrame
        tsv file as a pandas dataframe
    """
    cfg_opt = get_input_config('input')
    database_dir = cfg_opt['cosmic_dir']  # COSMIC_nuc database directory
    gene_dir = gene_name[0].upper() + '/'  # gene tsv in alphabetical directory listing
    tsv_path = database_dir + gene_dir + gene_name + '.tsv'  # path to tsv file
    df = pd.read_csv(tsv_path, sep='\t')
    return df


def drop_table(tbl_name,
               genes_db_path='',
               kind='sqlite'):
    """Drop a table from database if exists.

    **Note:** This function was written because pandas has a bug.
    If pandas was working then the write_frame method could just
    replace existing contents with out the need for me to drop the
    table. The bug is found here:
        https://github.com/pydata/pandas/issues/2971

    Parameters
    ----------
    tbl_name : str
        name of table to drop
    kind : str, ['sqlite' | 'mysql']
        type of database
    """
    if not genes_db_path:
        # if db not specified, use config file
        genes_db_path = get_db_config('2020plus')['db']
    if kind == 'sqlite':
        with sqlite3.connect(genes_db_path) as cur:
            sql = "DROP TABLE IF EXISTS %s" % tbl_name
            cur.execute(sql)
    elif kind == 'mysql':
        with get_cosmic_db() as cur:
            sql = "DROP TABLE IF EXISTS %s" % tbl_name
            cur.execute(sql)


def create_empty_table(tbl_name, db_path, colnames, coltypes):
    # drop table if exists
    drop_table(tbl_name, db_path, kind='sqlite')

    # make empty maf_mutation table
    conn = sqlite3.connect(db_path)  # open connection
    cur = conn.cursor()
    col_info_list = [' '.join(x) for x in zip(colnames, coltypes)]
    col_info_str = ', '.join(col_info_list)
    sql = "CREATE TABLE {0}({1});".format(tbl_name, col_info_str)
    cur.execute(sql)
    conn.commit()


def get_gene_length():
    # get db
    db_path = get_db_config('2020plus')['db']

    # query for gene length
    conn = sqlite3.connect(db_path)
    sql = "SELECT gene, gene_length FROM gene_features"
    df = psql.frame_query(sql, con=conn)
    df = df.set_index('gene')
    conn.close()

    return df


def cummin(x):
    """A python implementation of the cummin function in R"""
    for i in range(1, len(x)):
        if x[i-1] < x[i]:
            x[i] = x[i-1]
    return x


def bh_fdr(pval):
    """A python implementation of the Benjamani-Hochberg FDR method.

    This code should always give precisely the same answer as using
    p.adjust(pval, method="BH") in R.

    Parameters
    ----------
    pval : list or array
        list/array of p-values

    Returns
    -------
    pval_adj : np.array
        adjusted p-values according the benjamani-hochberg method
    """
    pval_array = np.array(pval)
    sorted_order = np.argsort(pval_array)
    original_order = np.argsort(sorted_order)
    pval_array = pval_array[sorted_order]

    # calculate the needed alpha
    n = float(len(pval))
    pval_adj = np.zeros(n)
    i = np.arange(1, n+1, dtype=float)[::-1]  # largest to smallest
    pval_adj = np.minimum(1, cummin(n/i * pval_array[::-1]))[::-1]
    return pval_adj[original_order]


# set up vogelstein oncogenes/tsgs
oncogene_list = read_oncogenes()
tsg_list = read_tsgs()
oncogene_set = set(oncogene_list)
tsg_set = set(tsg_list)

# significantly mutate genes from kandoth et al
#smg_list = read_smg()
#smg_set = set(smg_list)

# cancer gene census
#cgc_list = read_cgc()
#cgc_set = set(cgc_list)

# olfactory receptors from mutsigcv
olfactory_list = read_olfactory_receptors()
olfactory_set = set(olfactory_list)
