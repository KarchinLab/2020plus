import utils.python
from utils.python.cosmic_db import get_cosmic_db
from utils.python.amino_acid import AminoAcid
import dna_substitutions
import utils.python.util as _utils
import plot_data
import mutation_types
import single_gene
import tables.sample_name as sample_name
import tables.cosmic_aa as cosmic_aa
import tables.cosmic_genomic as cosmic_genomic
import missense
import pandas.io.sql as psql
import sqlite3
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

    df = psql.frame_query("SELECT * FROM nucleotide", con=conn)  # get all
    mtypes = _utils.get_mutation_types(df['AminoAcid'])
    df['mut_types'] = mtypes  # add mutation types to SQL output
    gene_to_indexes = df.groupby('Gene').groups

    # aggregate info
    design_matrix = []
    not_used_types = ['not valid',
                      'missing',
                      'unknown effect']  # only include known mutations
    for gene, indexes in gene_to_indexes.iteritems():
        tmp_df = df.ix[indexes]
        gene_pos_counter = {}
        mut_type_ctr = OrderedDict([['missense', 0],
                                    ['frame shift', 0],
                                    ['synonymous', 0],
                                    #['not valid', 0],
                                    ['indel', 0],
                                    #['missing', 0],
                                    ['nonsense', 0]])
                                    #['unknown effect', 0]])
        for hgvs in tmp_df['AminoAcid']:
            aa = AminoAcid(hgvs)
            if aa.mutation_type not in not_used_types:
                # do not use 'missing', 'unkown effect' or 'not valid'
                if aa.mutation_type == 'missense':
                    # keep track of missense pos for recurrency
                    gene_pos_counter.setdefault(aa.pos, 0)
                    gene_pos_counter[aa.pos] += 1
                mut_type_ctr[aa.mutation_type] += 1

        # needs to have at least one count
        if sum(mut_type_ctr.values()):
            recurrent_cts = sum([cts for cts in gene_pos_counter.values() if cts > 1])
            mut_type_ctr['missense'] -= recurrent_cts  # subtract off the recurrent missense
            design_matrix.append([gene, recurrent_cts] + list(mut_type_ctr.values()))
    header = [['gene', 'recurrent missense'] + list(mut_type_ctr)]
    logger.info('Finished creating design matrix.')
    return header + design_matrix


def count_gene_mutations(conn):
    logger.info('Counting number of mutations for each gene . . .')

    sql = """SELECT Gene, COUNT(*) as count
          FROM nucleotide GROUP BY Gene
          ORDER BY count DESC;"""
    logger.debug('Gene mutation count SQL statement: ' + sql)

    df = psql.frame_query(sql, con=conn)  # execute query
    logger.info('Finished getting gene mutation counts.')
    return df


def cosmic_nuc_main():
    cfg_opts = _utils.get_output_config('stats')  # get config
    conn = get_cosmic_db()

    # check info about COSMIC_nuc tables
    cosmic_aa.main()  # check cosmic_aa table
    sample_name.main()  # info related to mutations for each sample
    cosmic_genomic.main()  # check cosmic_genomic table

    # get design matrix
    design_matrix = generate_design_matrix(conn)
    with open(_utils.result_dir + cfg_opts['gene_design_matrix'], 'wb') as handle:
        csv.writer(handle, delimiter='\t').writerows(design_matrix)
    plot_data.pca_plot(_utils.result_dir + cfg_opts['gene_design_matrix'],
                       _utils.plot_dir + cfg_opts['pca_plot'])

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

    # handle protein missense mutations
    missense.main()

    # look at individual genes
    with open(_utils.config_dir + 'single_gene.txt') as handle:
        for row in handle:
            gene = row.strip()
            try:
                single_gene.main(gene)
            except:
                # be careful this catches all exceptions which might
                # hide other errors
                logger.debug('(Problem) Gene not found: %s' % gene)

    conn.close()  # close COSMIC_nuc connection


def main(db):
    cfg_opts = _utils.get_output_config('stats')  # get config

    if db == 'cosmic_nuc':
        conn = get_cosmic_db()  # connect to COSMIC_nuc

        # check info about COSMIC_nuc tables
        cosmic_aa.main()  # check cosmic_aa table
        cosmic_genomic.main()  # check cosmic_genomic table
    else:
        # connect to sqlite db at data/genes.db
        genes_db_path = _utils.get_db_config('genes')['db']
        conn = sqlite3.connect(genes_db_path)

    sample_name.main(conn)  # info related to mutations for each sample

    # get design matrix
    design_matrix = generate_design_matrix(conn)
    with open(_utils.result_dir + cfg_opts['gene_design_matrix'], 'wb') as handle:
        csv.writer(handle, delimiter='\t').writerows(design_matrix)
    plot_data.pca_plot(_utils.result_dir + cfg_opts['gene_design_matrix'],
                       _utils.plot_dir + cfg_opts['pca_plot'])

    # handle DNA substitutions
    dna_substitutions.main(conn)

    # handle stats related to mutation types
    mutation_types.main(conn)

    # gene mutation counts
    gene_ct_df = count_gene_mutations(conn)
    gene_ct_df.set_index('Gene').to_csv(_utils.result_dir + 'gene_mutation_counts.txt',
                                        sep='\t')
    plot_data.gene_mutation_histogram(gene_ct_df['count'],
                                      _utils.plot_dir + 'gene_mutations.histogram.png')
    plot_data.cumulative_gene_mutation(gene_ct_df['count'],
                                       _utils.plot_dir + 'gene_mutation.cumulative.png')

    # look at individual genes
    with open(_utils.config_dir + 'single_gene.txt') as handle:
        for row in handle:
            gene = row.strip()
            try:
                single_gene.main(gene, conn)
            except:
                # be careful this catches all exceptions which might
                # hide other errors
                logger.debug('(Problem) Gene not found: %s' % gene)

    conn.close()  # close data/genes.db connection
