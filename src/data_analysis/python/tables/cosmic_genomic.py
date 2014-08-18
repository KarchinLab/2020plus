"""
Checks the cosmic_genomic table in the COSMIC_nuc databse
particularly for mutation types present in the table.
"""

import pandas.io.sql as psql
import utils.python.util as _utils
from utils.python.cosmic_db import get_cosmic_db
import data_analysis.python.plot_data as plot_data
import logging

logger = logging.getLogger(__name__)

def count_amino_acids(conn):
    """Count the amino acid mutation types (missense, indel, etc.).

    NOTE: this function only counts the number of rows for a mutation
    type in the cosmic_genomic table. It does not include occurence
    information from the occurrence column.
    """
    df = psql.frame_query("SELECT aachange FROM `cosmic_genomic`", con=conn)
    unique_cts = _utils.count_mutation_types(df['aachange'])
    return unique_cts


def count_nucleotides(conn):
    """Count the number of nucleotide substitutions.
    """
    sql = ("SELECT refbase as initial, altbase as mutated, SUM(occurrences) as count "
          "FROM cosmic_genomic "
          "GROUP BY refbase, altbase "
          "ORDER BY count DESC")
    df = psql.frame_query(sql, con=conn)
    return df


def main():
    cfg_opts = _utils.get_output_config('cosmic_genomic')
    out_dir = _utils.result_dir  # output dir for txt files
    plot_dir = _utils.plot_dir  # output dir for figures

    conn = get_cosmic_db()  # connect to COSMIC_nuc

    # examine mutation types for the cosmic_aa table
    mut_cts = count_amino_acids(conn)  # all mutation cts
    mut_cts.to_csv(out_dir + cfg_opts['aa_type'], sep='\t')
    plot_data.mutation_types_barplot(mut_cts,
                                     save_path=plot_dir + cfg_opts['aa_type_barplot'],
                                     title=r'Protein Mutations by Type for '
                                           r'\textit{cosmic\_genomic} Table')

    # count dna substitutions
    nuc_counts = count_nucleotides(conn)
    nuc_save_path = out_dir + cfg_opts['sub_out']
    nuc_counts.to_csv(nuc_save_path, sep='\t')

    # plot dna substitutions
    all_heatmap = _utils.plot_dir + cfg_opts['sub_heatmap']
    plot_data.nuc_substitution_heatmap(nuc_save_path, all_heatmap)
    all_barplot = _utils.plot_dir + cfg_opts['sub_barplot']
    plot_data.nuc_substitution_barplot(nuc_save_path, all_barplot,
                                       title=r'DNA Substitution Mutations in '
                                             r'\textit{cosmic\_genomic} Table')
    conn.close()  # close connection to COSMIC_nuc
