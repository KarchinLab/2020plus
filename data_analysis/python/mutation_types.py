"""
The mutation_types module stratifies counts by mutation types
for amino acids (missense, indel, frame shift, nonsense, and synonymous)
and nucleotides (substitution, insertions, and deletions).
"""

import utils.python.util as _utils
import pandas as pd
import pandas.io.sql as psql
import logging


def count_amino_acids(conn):
    """Count the amino acid mutation types (missense, indel, etc.).
    """
    df = psql.frame_query("""SELECT * FROM `nucleotide`""", con=conn)
    unique_cts = _utils.count_mutation_types(df['AminoAcid'])
    return unique_cts


def count_nucleotides(conn):
    """Count the nucleotide mutation types (substitution, indels)
    """
    sql = "SELECT Nucleotide FROM `nucleotide`"
    df = psql.frame_query(sql, con=conn)
    unique_cts = _utils.count_mutation_types(df['Nucleotide'], kind='nucleotide')
    return unique_cts


def count_oncogenes(conn):
    logger = logging.getLogger(__name__)
    logger.info('Counting oncogene mutation types . . .')

    # prepare sql statement
    oncogenes = _utils.read_oncogenes()
    sql = "SELECT * FROM `nucleotide` WHERE Gene in " + str(oncogenes)
    logger.debug('Oncogene SQL statement: ' + sql)

    df = psql.frame_query(sql, con=conn)  # execute query

    # count mutation types
    aa_counts = _utils.count_mutation_types(df['AminoAcid'])
    nuc_counts = _utils.count_mutation_types(df['Nucleotide'],
                                             kind='nucleotide')
    logger.info('Finished counting oncogene mutation types.')
    return aa_counts, nuc_counts


def count_tsg(conn):
    logger = logging.getLogger(__name__)
    logger.info('Counting tumor suppressor gene mutation types . . .')

    # prepare sql statement
    tsgs = _utils.read_tsgs()
    sql = "SELECT * FROM `nucleotide` WHERE Gene in " + str(tsgs)
    logger.debug('Oncogene SQL statement: ' + sql)

    df = psql.frame_query(sql, con=conn)  # execute query

    # count mutation types
    aa_counts = _utils.count_mutation_types(df['AminoAcid'])
    nuc_counts = _utils.count_mutation_types(df['Nucleotide'],
                                             kind='nucleotide')
    logger.info('Finished counting tumor suppressor gene mutation types.')
    return aa_counts, nuc_counts


def count_gene_types(file_path='data_analysis/results/gene_design_matrix.txt'):
    """Returns protein mutation type counts by gene type (oncogenes, tsg, other).

    Kwargs:
        file_path (str): path to mutation type cts by gene file

    Returns:
        pd.DataFrame: mutation type counts by gene type
    """
    df = pd.read_csv(file_path, sep='\t')
    df['gene_type'] = df['gene'].apply(_utils.classify_gene)
    mut_ct_df = df.iloc[:, 1:]  # remove the "gene" column
    mut_ct_df = mut_ct_df.groupby('gene_type').sum()  # get counts for each gene type
    return mut_ct_df
