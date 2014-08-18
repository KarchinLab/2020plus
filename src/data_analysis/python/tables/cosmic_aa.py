"""
The cosmic_aa module checks the characteristics of the
cosmic_aa table in the COSMIC_nuc database. For example,
how many of the amino acid changes are missense? This
could potentially show whether other information was
filtered from the table.
"""

import pandas.io.sql as psql
import src.utils.python.util as _utils
from src.utils.python.cosmic_db import get_cosmic_db
import src.data_analysis.python.plot_data as plot_data
import logging

logger = logging.getLogger(__name__)

def count_amino_acids(conn):
    """Count the amino acid mutation types (missense, indel, etc.).

    NOTE: this function only counts the number of rows for a mutation
    type in the cosmic_genomic table. It does not include occurence
    information from the occurrence column.
    """
    df = psql.frame_query("SELECT aachange FROM `cosmic_aa`", con=conn)
    unique_cts = _utils.count_mutation_types(df['aachange'])
    return unique_cts


def main():
    cfg_opts = _utils.get_output_config('cosmic_aa')
    out_dir = _utils.result_dir  # output dir for txt files
    plot_dir = _utils.plot_dir  # output dir for figures

    conn = get_cosmic_db()  # connect to COSMIC_nuc

    # examine mutation types for the cosmic_aa table
    mut_cts = count_amino_acids(conn)  # all mutation cts
    mut_cts.to_csv(out_dir + cfg_opts['aa_type'], sep='\t')
    plot_data.mutation_types_barplot(mut_cts,
                                     save_path=plot_dir + cfg_opts['aa_type_barplot'],
                                     title=r'Protein Mutations by Type for '
                                           r'\textit{cosmic\_aa} Table')

    conn.close()  # close connection to COSMIC_nuc
