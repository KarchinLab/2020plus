import utils.python
from utils.python.cosmic_db import get_cosmic_db
from utils.python.amino_acid import AminoAcid
import dna_substitutions
import utils.python.util as _utils
import plot_data
import mutation_types
import sample_name
import missense
import tp53
import pandas.io.sql as psql
import csv
from collections import OrderedDict
import logging

logger = logging.getLogger(__name__)

def generate_design_matrix(conn):
    """Generate a design matrix potentially useful for classifying genes.

    Args:
        conn (MySQLdb connection): database connection to mysql

    Returns:
        pd.DataFrame: Design matrix
    """
    logger.info('Creating design matrix . . .')

    df = psql.frame_query("SELECT * FROM `nucleotide`", con=conn)  # get all
    mtypes = _utils.get_mutation_types(df['AminoAcid'])
    df['mut_types'] = mtypes  # add mutation types to SQL output
    gene_to_indexes = df.groupby('Gene').groups

    # aggregate info
    design_matrix = []
    for gene, indexes in gene_to_indexes.iteritems():
        tmp_df = df.ix[indexes]
        gene_pos_counter = {}
        mut_type_ctr = OrderedDict([['missense', 0],
                                    ['frame shift', 0],
                                    ['synonymous', 0],
                                    ['not valid', 0],
                                    ['indel', 0],
                                    ['missing', 0],
                                    ['nonsense', 0]])
        for hgvs in tmp_df['AminoAcid']:
            aa = AminoAcid(hgvs)
            if aa.mutation_type == 'missense':
                gene_pos_counter.setdefault(aa.pos, 0)
                gene_pos_counter[aa.pos] += 1
            mut_type_ctr[aa.mutation_type] += 1

        recurrent_cts = sum([cts for cts in gene_pos_counter.values() if cts > 1])
        mut_type_ctr['missense'] -= recurrent_cts  # subtract off the recurrent missense
        design_matrix.append([gene, recurrent_cts] + list(mut_type_ctr.values()))
    header = [['gene', 'recurrent missense'] + list(mut_type_ctr)]
    return header + design_matrix


def count_gene_mutations(conn):
    logger.info('Counting number of mutations for each gene . . .')

    sql = """SELECT Gene, COUNT(*) as count
          FROM `nucleotide` GROUP BY Gene
          ORDER BY count DESC;"""
    logger.debug('Gene mutation count SQL statement: ' + sql)

    df = psql.frame_query(sql, con=conn)  # execute query
    logger.info('Finished getting gene mutation counts.')
    return df


def main():
    # count mutation types
    conn = get_cosmic_db()

    # look at TP53
    tp53.main()

    """
    # get design matrix
    design_matrix = generate_design_matrix(conn)
    with open(_utils.result_dir + 'gene_design_matrix.txt', 'wb') as handle:
        csv.writer(handle, delimiter='\t').writerows(design_matrix)
    plot_data.pca_plot()

    # handle DNA substitutions
    dna_substitutions.main()

    # handle stats related to mutation types
    mutation_types.main()

    # gene mutation counts
    gene_ct_df = count_gene_mutations(conn)
    gene_ct_df.set_index('Gene').to_csv(_utils.result_dir + 'gene_mutation_counts.txt',
                                        sep='\t')
    plot_data.gene_mutation_histogram(gene_ct_df['count'])
    plot_data.cumulative_gene_mutation(gene_ct_df['count'])

    # get information related to mutations for each sample
    sample_name.main()

    # handle protein missense mutations
    missense.main()
    """
    conn.close()


if __name__=="__main__":
    main()
