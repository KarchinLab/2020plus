"""
The sample_name module analyzes characteristics of samples in
the `COSMIC_nuc` database by using the `nucleotide`.SampleName
column.
"""

import pandas as pd
import pandas.io.sql as psql


def count_mutated_genes(conn):
    """Count the number of genes that are mutated in each sample.

    Args:
        conn (MySQLdb connection): connection to COSMIC_nuc

    Returns:
        pd.DataFrame: two column data frame of sample names and gene cts
    """
    sql = ('SELECT x.SampleName, SUM(x.gene_indicator) as GeneCounts'
          ' FROM ('
          '     SELECT SampleName, Gene, 1 as gene_indicator'
          '     FROM `nucleotide`'
          '     GROUP BY SampleName, Gene'
          ' ) x GROUP BY SampleName'
          ' ORDER BY GeneCounts Desc;')
    df = psql.frame_query(sql, con=conn)
    df.GeneCounts = df.GeneCounts.astype(int)  # pandas is not auto detecting
    return df


def count_mutations(conn):
    """Count the number of mutations in each sample.

    Args:
        conn (MySQLdb connection): connection to COSMIC_nuc

    Returns:
        pd.DataFrame: two column data frame of sample names and mutation cts
    """
    sql = ('SELECT x.SampleName, SUM(x.mut_indicator) as MutationCounts'
          ' FROM ('
          '     SELECT SampleName, 1 as mut_indicator'
          '     FROM `nucleotide`'
          ' ) x GROUP BY SampleName'
          ' ORDER BY MutationCounts Desc;')
    df = psql.frame_query(sql, con=conn)
    df.MutationCounts = df.MutationCounts.astype(int)  # pandas is not auto detecting
    return df
