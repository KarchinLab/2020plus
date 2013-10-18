import pandas as pd
from amino_acid import AminoAcid
import logging


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
    with open('data_analysis/gene_lists/oncogenes_vogelstein.txt', 'r') as handle:
        oncogenes = tuple(gene.strip() for gene in handle.readlines())
    return oncogenes


def read_tsgs():
    """Reads in the tumor suppressor genes from vogelsteins' science paper.

    Oncogenes from supplementary 2A:
    http://www.sciencemag.org/content/339/6127/1546.full

    Returns:
        tsgs (tuple): tuple of gene names considered as tumor suppressors
    """
    with open('data_analysis/gene_lists/tsg_vogelstein.txt', 'r') as handle:
        tsgs = tuple(gene.strip() for gene in handle.readlines())
    return tsgs


def get_mutation_types(hgvs_iterable):
    """Classify each protein HGVS mutation as a certain type.

    Args:
        hgvs_iterable (iterable): iterable container with HGVS mutaiton strings

    Returns:
        pd.Series: container of protein mutation types in same order as input
    """
    mut_type = []
    for hgvs_aa in hgvs_iterable:
        aa = AminoAcid(hgvs=hgvs_aa)
        mut_type.append(aa.mutation_type)
    mut_type_series = pd.Series(mut_type)
    return mut_type_series


def count_mutation_types(hgvs_iterable):
    """Count protein mutation types from HGVS strings (missense, indels, etc.).

    Args:
        hgvs_iterable (iterable): An iterable object containing protein HGVS

    Returns:
        pd.Series: A pandas series object counting protein mutation types
    """
    mut_type_series = get_mutation_types(hgvs_iterable)  # get mutation types
    unique_cts = mut_type_series.value_counts() # count mutation types
    return unique_cts

