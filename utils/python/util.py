import pandas as pd
from amino_acid import AminoAcid
from nucleotide import Nucleotide
import ConfigParser
import logging

# setup directory paths
plot_dir = 'data_analysis/plots/genes/'
result_dir = 'data_analysis/results/genes/'
clf_plot_dir = 'classify/plots/'
clf_result_dir = 'classify/results/'
config_dir = 'config/'

def read_aa_properties(file_path):
    """Read aa property counts from the data_analysis/results folder.

    Args:
        file_path (str): path to aa_change.properties.txt

    Returns:
        DataFrame. contains mutation counts for amino acid chemical properties
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

    Returns:
        oncogenes (tuple): tuple of gene names considered oncogenes
    """
    with open('data/gene_lists/oncogenes_vogelstein.txt', 'r') as handle:
        oncogenes = tuple(gene.strip() for gene in handle.readlines())
    return oncogenes


def read_tsgs():
    """Reads in the tumor suppressor genes from vogelsteins' science paper.

    Oncogenes from supplementary 2A:
    http://www.sciencemag.org/content/339/6127/1546.full

    Returns:
        tsgs (tuple): tuple of gene names considered as tumor suppressors
    """
    with open('data/gene_lists/tsg_vogelstein.txt', 'r') as handle:
        tsgs = tuple(gene.strip() for gene in handle.readlines())
    return tsgs


def read_smg():
    """Reads in the significantly mutated genes from kandoth et al's
    nature paper.

    The paper was from the Pan-Cancer effort from TCGA.

    SMGs from supplementary:
    http://www.nature.com/nature/journal/v502/n7471/full/nature12634.html

    Returns:
        smgs (tuple): tuple of gene names considered as significantly mutated
    """
    with open('data/gene_lists/smg_kandoth.txt', 'r') as handle:
        smgs = tuple(gene.strip() for gene in handle.readlines())
    return smgs


def classify_gene(gene):
    """Return whether the gene is an oncogene, tsg, or other.

    Args:
        gene (str): Official gene name

    Returns:
        Str: 'oncogene', 'tsg', or 'other'
    """
    if gene in oncogene_set:
        return 'oncogene'
    elif gene in tsg_set:
        return 'tsg'
    else:
        return 'other'


def get_mutation_types(hgvs_iterable, kind='amino acid'):
    """Classify each protein HGVS mutation as a certain type.

    Args:
        hgvs_iterable (iterable): iterable container with HGVS mutaiton strings

    Returns:
        pd.Series: container of protein mutation types in same order as input
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

    Args:
        hgvs_iterable (iterable): An iterable object containing protein HGVS

    Returns:
        pd.Series: A pandas series object counting protein mutation types
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


def get_input_config(section):
    """Returns the config object to input.cfg."""
    cfg = ConfigParser.ConfigParser()
    cfg.read(config_dir + 'input.cfg')
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

    Args:
        gene_name (str): gene name

    Returns:
        pd.DataFrame: tsv file as a pandas dataframe
    """
    cfg_opt = get_input_config('input')
    database_dir = cfg_opt['cosmic_dir']  # COSMIC_nuc database directory
    gene_dir = gene_name[0].upper() + '/'  # gene tsv in alphabetical directory listing
    tsv_path = database_dir + gene_dir + gene_name + '.tsv'  # path to tsv file
    df = pd.read_csv(tsv_path, sep='\t')
    return df


# set up vogelstein oncogenes/tsgs
oncogene_list = read_oncogenes()
tsg_list = read_tsgs()
oncogene_set = set(oncogene_list)
tsg_set = set(tsg_list)

# significantly mutate genes from kandoth et al
smg_list = read_smg()
smg_set = set(smg_list)
