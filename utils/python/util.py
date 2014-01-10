import pandas as pd
from amino_acid import AminoAcid
from nucleotide import Nucleotide
from cosmic_db import get_cosmic_db
import sqlite3
import pandas.io.sql as psql
import ConfigParser
import logging

logger = logging.getLogger(__name__)

config_dir = 'config/'

def get_input_config(section):
    """Returns the config object to input.cfg."""
    cfg = ConfigParser.ConfigParser()
    cfg.read(config_dir + 'input.cfg')
    cfg_options = dict(cfg.items(section))
    return cfg_options

# setup directory paths
_opts = get_input_config('result')
plot_dir = _opts['plot_dir']
result_dir = _opts['result_dir']
clf_plot_dir = _opts['clf_plot_dir']
clf_result_dir = _opts['clf_result_dir']

def read_aa_properties(file_path):
    """Read aa property counts from the data_analysis/results folder.

    **Parameters**

    file_path : str
        path to aa_change.properties.txt

    **Returns**

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

    **Returns**

    oncogenes : tuple
        tuple of gene names considered oncogenes
    """
    cfg_opts = get_input_config('input')
    with open(cfg_opts['oncogene'], 'r') as handle:
        oncogenes = tuple(gene.strip() for gene in handle.readlines())
    return oncogenes


def read_tsgs():
    """Reads in the tumor suppressor genes from vogelsteins' science paper.

    Oncogenes from supplementary 2A:
    http://www.sciencemag.org/content/339/6127/1546.full

    **Returns**

    tsgs : tuple
        tuple of gene names considered as tumor suppressors
    """
    cfg_opts = get_input_config('input')
    with open(cfg_opts['tsg'], 'r') as handle:
        tsgs = tuple(gene.strip() for gene in handle.readlines())
    return tsgs


def read_smg():
    """Reads in the significantly mutated genes from kandoth et al's
    nature paper.

    The paper was from the Pan-Cancer effort from TCGA.

    SMGs from supplementary:
    http://www.nature.com/nature/journal/v502/n7471/full/nature12634.html

    **Returns**

    smgs : tuple
        tuple of gene names considered as significantly mutated
    """
    cfg_opts = get_input_config('input')
    with open(cfg_opts['smg'], 'r') as handle:
        smgs = tuple(gene.strip() for gene in handle.readlines())

    # open connection
    try:
        # if DB is not created this will throw an error
        gene_db_path = get_db_config('genes')['db']
        conn = sqlite3.connect(gene_db_path)

        sql = ("SELECT DISTINCT Gene"
              " FROM nucleotide"
              " WHERE Gene in " + str(smgs))
        df = psql.frame_query(sql, con=conn)
        conn.close()  # close connection

        # get significantly mutated genes found in database
        smgs_in_database = tuple(df['Gene'])
        logger.info('There are only %d/%d significantly mutated genes found in the database.'
                    % (len(smgs_in_database), len(smgs)))
    except:
        smgs_in_database = smgs
    return smgs_in_database


def read_olfactory_receptors():
    """Reads in the significant olfactory receptors from Mutsigcv.

    **Returns**

    olfactory : tuple
        tuple of gene names considered as olfactory receptors
    """
    cfg_opts = get_input_config('input')
    with open(cfg_opts['olfactory_receptors'], 'r') as handle:
        olfactory = tuple(gene.strip() for gene in handle.readlines())

    # open connection
    try:
        # it table is not found just catch exception
        gene_db_path = get_db_config('genes')['db']
        conn = sqlite3.connect(gene_db_path)

        sql = ("SELECT DISTINCT Gene"
              " FROM nucleotide"
              " WHERE Gene in " + str(olfactory))

        df = psql.frame_query(sql, con=conn)
        conn.close()  # close connection

        # get significantly mutated genes found in database
        olfactory_in_database = tuple(df['Gene'])
        logger.info('There are only %d/%d olfactory receptors found in the database.'
                    % (len(olfactory_in_database), len(olfactory)))
    except:
        olfactory_in_database = olfactory
    return olfactory_in_database


def classify_gene(gene):
    """Return whether the gene is an oncogene, tsg, or other.

    **Parameters**

    gene : str
        Official gene name

    **Returns**

        Str, ['oncogene' | 'tsg' | 'other']
    """
    if gene in oncogene_set:
        return 'oncogene'
    elif gene in tsg_set:
        return 'tsg'
    else:
        return 'other'


def get_mutation_types(hgvs_iterable, kind='amino acid'):
    """Classify each protein HGVS mutation as a certain type.

    **Parameters**

    hgvs_iterable : iterable
        iterable container with HGVS mutaiton strings

    **Returns**

    mut_type_series : pd.Series
        container of protein mutation types in same order as input
    """
    mut_type = []
    if kind == 'amino acid':
        for hgvs_aa in hgvs_iterable:
            aa = AminoAcid(hgvs=hgvs_aa)
            mut_type.append(aa.mutation_type)
    elif kind == 'nucleotide':
        for hgvs_nuc in hgvs_iterable:
            nuc = Nucleotide(hgvs=hgvs_nuc)
            mut_type.append(nuc.mutation_type)
    mut_type_series = pd.Series(mut_type)
    return mut_type_series


def count_mutation_types(hgvs_iterable, kind='amino acid'):
    """Count mutation types from HGVS protein strings (missense, indels, etc.)
    and DNA strings (substitutions, indels).

    **Parameters**

    hgvs_iterable : iterable
        An iterable object containing protein HGVS

    **Returns**

    unique_cts : pd.Series
        A pandas series object counting protein mutation types
    """
    mut_type_series = get_mutation_types(hgvs_iterable, kind=kind)  # get mutation types
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

    **Parameters**

    gene_name : str
        gene name

    **Returns**

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
               kind='sqlite'):
    """Drop a table from database if exists.

    **Note:** This function was written because pandas has a bug.
    If pandas was working then the write_frame method could just
    replace existing contents with out the need for me to drop the
    table. The bug is found here:
        https://github.com/pydata/pandas/issues/2971

    **Parameters**

    tbl_name : str
        name of table to drop
    kind : str, ['sqlite' | 'mysql']
        type of database
    """
    genes_db_path = get_db_config('genes')['db']
    if kind == 'sqlite':
        with sqlite3.connect(genes_db_path) as cur:
            sql = "DROP TABLE IF EXISTS %s" % tbl_name
            cur.execute(sql)
    elif kind == 'mysql':
        with get_cosmic_db() as cur:
            sql = "DROP TABLE IF EXISTS %s" % tbl_name
            cur.execute(sql)


# set up vogelstein oncogenes/tsgs
oncogene_list = read_oncogenes()
tsg_list = read_tsgs()
oncogene_set = set(oncogene_list)
tsg_set = set(tsg_list)

# significantly mutate genes from kandoth et al
smg_list = read_smg()
smg_set = set(smg_list)

# olfactory receptors from mutsigcv
olfactory_list = read_olfactory_receptors()
olfactory_set = set(olfactory_list)
