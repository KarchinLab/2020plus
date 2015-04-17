"""The gene_features module generates a DB table for features
like gene length and misc. features from the MutSigCV paper (fig s5).

The MutSigCV paper can be found here:
    http://www.nature.com/nature/journal/v499/n7457/full/nature12213.html

The MutSigCV paper suggested that the background mutation rate for genes is important
for identifying statistically significant cancer genes. However, it is
not clear how important those features are for supervised learning on
"known" oncogenes and tsgs.
"""

import src.utils.python.util as _utils
import pandas as pd
import pandas.io.sql as psql
import sqlite3
import string
import os
import logging

logger = logging.getLogger(__name__)

def cnv_gain(db_path):
    logger.info('Calculating CNV gain features . . .')
    conn = sqlite3.connect(db_path)  # open connection
    sql = ('SELECT x.gene as gene, SUM(x.cnv_indicator) as cnv_gain '
           'FROM( '
           '    SELECT gf.gene, 1 as cnv_indicator'
           '    FROM cosmic_cnv as cc, gene_features as gf'
           '    WHERE cc.MUT_TYPE="gain" AND (cc.Chromosome=gf.chr AND (cc.G_Start<=gf.start AND cc.G_Stop>=gf.end))'
           ') x GROUP BY x.gene')
    cnv_gain_df = psql.frame_query(sql, con=conn)
    conn.close()
    logger.info('Finished calculating CNV gain features.')
    return cnv_gain_df


def cnv_loss(db_path):
    logger.info('Calculating CNV loss features . . .')
    conn = sqlite3.connect(db_path)  # open connection
    sql = ('SELECT x.gene as gene, SUM(x.cnv_indicator) as cnv_loss '
           ' FROM( '
           '     SELECT gf.gene, 1 as cnv_indicator '
           '     FROM cosmic_cnv as cc, gene_features as gf '
           '     WHERE cc.MUT_TYPE="loss" AND (cc.Chromosome=gf.chr AND ((cc.G_Start<=gf.start AND cc.G_Stop>=gf.end) OR (cc.G_Start<gf.end AND cc.G_Start>gf.start) OR (cc.G_Stop<gf.end AND cc.G_STOP>gf.start))) '
           ' ) x GROUP BY x.gene')
    cnv_loss_df = psql.frame_query(sql, con=conn)
    conn.close()
    logger.info('Finished calculating CNV loss features.')
    return cnv_loss_df


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
        print mydir
        dir_path = fasta_dir + mydir + '/'
        for file_name in os.listdir(dir_path):
            if '_protein' in file_name and '_ENST' not in file_name:
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
    db_opts = _utils.get_db_config('2020plus')

    # get absolute path for cosmic data
    cosmic_path = os.path.join(_utils.proj_dir, in_opts['cosmic_path'])

    # get data for gene_features table
    logger.info('Processing features for gene_features table ...')
    if os.path.isdir(cosmic_path):
        gene_length = recursive_gene_length(in_opts['fasta_dir'])
        genes, lengths = zip(*gene_length.items())
        gene_length_df = pd.DataFrame({'gene': genes, 'gene length': lengths})
    else:
        gene_length_df = pd.read_csv(cosmic_path, sep='\t')
        gene_length_df = gene_length_df[['Gene name', 'Gene CDS length']]
        gene_length_df = gene_length_df.rename(columns={'Gene name': 'gene',
                                                        'Gene CDS length': 'gene length'})
        gene_length_df.drop_duplicates(cols=['gene'], inplace=True)

    # merge in data from mutsig and biogrid
    mutsigcv_feature_path = os.path.join(_utils.proj_dir, in_opts['mutsigcv_features'])
    df = pd.read_csv(mutsigcv_feature_path, sep='\t')
    df = pd.merge(gene_length_df, df, how='left', on='gene')  # merge the data frames
    biogrid_path = os.path.join(_utils.proj_dir, 'data/biogrid_stats.txt')
    biogrid_df = pd.read_csv(biogrid_path, sep='\t')
    df = pd.merge(df, biogrid_df, how='left', on='gene')

    # path to database
    db_path = db_path if db_path else db_opts['db']

    # get info on cnv's overlapping genes
    #tmp_gain = cnv_gain(db_path)
    #df = pd.merge(df, tmp_gain, on='gene', how='left')
    #df['cnv_gain'] = df['cnv_gain'].fillna(0)
    #tmp_loss = cnv_loss(db_path)
    #df = pd.merge(df, tmp_loss, on='gene', how='left')
    #df['cnv_loss'] = df['cnv_loss'].fillna(0)

    #df['gene length'] = gene_length_series
    #df['gene'] = df.index  # add gene names as a column (not just an index)

    logger.info('Finished processing features for gene_features table.')

    # save database
    save_db(df, db_path)
