"""The gene_features module generates a DB table for features
like gene length and misc. features from the MutSigCV paper.

The MutSigCV paper can be found here:
    http://www.nature.com/nature/journal/v499/n7457/full/nature12213.html

All features except gene length are taken from Table S5. The MutSigCV
paper suggested that the background mutation rate for genes is important
for identifying statistically significant cancer genes. However, it is
not clear how important those features are for supervised learning on
"known" oncogenes and tsgs.
"""

import utils.python.util as _utils
import pandas as pd
import pandas.io.sql as psql
import sqlite3
import string
import os
import logging

logger = logging.getLogger(__name__)

def calc_gene_length(file_path):
    """Read in a FASTA file and calculate sequence length.

    Assumes a typical one line header for a FASTA file.

    **Parameters**

    file_path : str
        Path to FASTA file

    **Returns**

    seq_len : int
        length of gene
    """
    with open(file_path) as handle:
        handle.readline()  # skip FASTA header
        seq = handle.read()  # read file into single string
        seq = seq.replace('\n', '')  # remove line breaks
    seq_len = len(seq)
    return seq_len


def recursive_gene_length(fasta_dir):
    """Recursively scans the FASTA directory to calc gene lengths.

    NOTE: assumes directories are ['0-9', 'A', .., 'Z']

    **Parameters**

    fasta_dir : str
        path to fasta directory downloaded from COSMIC

    **Returns**

    gene_length_dict : dict
        keys=gene name, values=gene length
    """
    logger.info('Recursively calculating length in FASTA directories . . .')
    gene_length_dict = {}
    mydirs = ['0-9'] + list(string.ascii_uppercase)
    for mydir in mydirs:
        dir_path = fasta_dir + mydir + '/'
        for file_name in os.listdir(dir_path):
            if '_protein' in file_name:
                gene_name = file_name.strip('_protein.txt')
                gene_length = calc_gene_length(dir_path + file_name)
                gene_length_dict[gene_name] = gene_length
    logger.info('Finished counting gene length.')
    return gene_length_dict


def save_db(df, genedb_path):
    """Saves the data into the gene_features table.

    If the table already exists, the table is droped and then
    re-inserted.

    **Parameters**

    df : pd.DataFrame
        data to insert into DB table
    genedb_path : str
        path to sqlite db
    """
    logger.debug('Dropping gene_features table IF EXISTS.')
    _utils.drop_table('gene_features', genes_db_path=genedb_path, kind='sqlite')  # drop table if exists
    logger.debug('After dropping gene_features table IF EXISTS.')

    logger.info('Saving gene_features table ...')
    conn = sqlite3.connect(genedb_path)  # open connection
    # save to sqlite3 database
    psql.write_frame(df,  # pandas dataframe
                     'gene_features',  # table name
                     con=conn,  # connection
                     flavor='sqlite',  # use sqlite
                     if_exists='replace')  # drop table if exists
    conn.close()
    logger.info('Finished saving gene_features table.')


def main(db_path):
    # get config files
    in_opts = _utils.get_input_config('input')
    db_opts = _utils.get_db_config('genes')

    # get data for gene_features table
    logger.info('Processing features for gene_features table ...')
    gene_length = recursive_gene_length(in_opts['fasta_dir'])
    genes, lengths = zip(*gene_length.items())
    gene_length_df = pd.DataFrame({'gene': genes, 'gene length': lengths})
    df = pd.read_csv(in_opts['mutsigcv_features'], sep='\t')
    df = pd.merge(gene_length_df, df, how='left', on='gene')  # merge the data frames
    #df['gene length'] = gene_length_series
    #df['gene'] = df.index  # add gene names as a column (not just an index)
    logger.info('Finished processing features for gene_features table.')

    db_path = db_path if db_path else db_opts['db']

    # save database
    save_db(df, db_path)
