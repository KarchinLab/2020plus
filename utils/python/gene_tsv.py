"""The gene_tsv module aggregates all of the tab delimited files
describing gene mutations that form the COSMIC_nuc database. It
saves the resulting tab delimited file into the data directory.
It subsequently makes a database from the file.

 * only lines with either a nucleotide or amino acid mutation
   are saved.
 * this module in basically just a script
 * it will likely break if paths/conventions change
 * mutations for alternative isoforms "_ENST..." are not used

Rows are filtered if:
 * either the amino acid or nucleotide hgvs column is empty
 * both columns indicate unkown effect c.? and p.?
 * SomaticStatus has word 'unkown'
"""

import utils.python.util as _utils
import pandas as pd
import pandas.io.sql as psql
import sqlite3
import string
import os


def skip_header(file_handle, skip_rows=8):
    """Skips the first "skip_row" lines of the file.

    To skip rows the next() method is called repeatedly.
    This skips the non-table part of the tab delimited gene files.

    Args
        file_handle (file): file object of tab delim file

    Kwargs
        skip_rows (int): number of lines to skip
    """
    for i in range(skip_rows):
        file_handle.next()
    return file_handle


def concatenate_genes(out_path, cosmic_dir):
    """Saves all information in tab delimited files from genes into
    a single file in the data directory.
    """
    # concatenate gene tab delim files
    with open(out_path, 'wb') as mywriter:
        # iterate through 'A'...'Z' directories
        FIRST_FLAG = True  # flag for writing header
        for letter in string.ascii_uppercase:
            for file_name in os.listdir(cosmic_dir + letter):
                if file_name.endswith('.tsv') and '_ENST' not in file_name:
                    # only use tab delim files that are
                    # not alternative isoforms
                    with open(cosmic_dir + letter + "/" + file_name) as handle:
                        if FIRST_FLAG:
                            # headers match the `Nucleotide` table in COSMIC_nuc
                            header = ['Gene', 'SampleName', 'COSMICSampleID',
                                      'AminoAcid', 'Nucleotide', 'PrimaryTissue',
                                      'Tissuesubtype1', 'Tissuesubtype2', 'Histology',
                                      'Histologysubtype1', 'Histologysubtype2', 'PubmedID',
                                      'studies', 'MutationID', 'SomaticStatus',
                                      'SampleSource', 'Zygosity', 'hg18chrom',
                                      'hg18start', 'hg18end', 'hg19chrom',
                                      'hg19start', 'hg19end']
                            mywriter.write('\t'.join(header) + '\n')  # write header
                            FIRST_FLAG = False

                        gene_name = file_name[:-4]  # file name before ".tsv"
                        handle = skip_header(handle)  # skip beginning lines
                        # add data
                        for line in handle:
                            split_line = line.split('\t')
                            if split_line[2] and split_line[3]:
                                # if line designates a mutation
                                if split_line[2] != 'p.?' or split_line[3] != 'c.?':
                                    # not unknown effect for AA and amino acid
                                    if 'unknown' not in split_line[-9].lower():
                                        # do not include unknown somatic status mutations
                                        mywriter.write(gene_name + "\t" + line)  # write with gene name

        # iterate through genes in the numeric directory
        numeric_dir = cosmic_dir + '0-9/'
        for file_name in os.listdir(numeric_dir):
            if file_name.endswith('.tsv') and '_ENST' not in file_name:
                # only use tab delim files that are
                # not alternative isoforms
                with open(numeric_dir + file_name) as handle:
                    gene_name = file_name[:-4]  # file name before ".tsv"
                    handle = skip_header(handle)  # skip beginning lines
                    for line in handle:
                        split_line = line.split('\t')
                        if split_line[2] and split_line[3]:
                            # if line designates a mutation
                            if split_line[2] != 'p.?' or split_line[3] != 'c.?':
                                # not uknown effect for both AA and nucleotide
                                if 'unknown' not in split_line[-9].lower():
                                    # do not include unknown somatic status mutations
                                    mywriter.write(gene_name + "\t" + line)  # write with gene name


def save_db(tsv_path, genedb_path):
    """Saves tab delim gene mutation file to a sqlite3 db.

    NOTE: Uses pandas to store all contents in memory and then
    saves to sqlite db. This may cause large memory usage.
    """
    df = pd.read_csv(tsv_path, sep='\t')  # read data

    # fix types that pandas gets wrong
    # see http://pandas.pydata.org/pandas-docs/dev/gotchas.html
    # for details on missing NA support for integers
    df['hg18chrom'] = df['hg18chrom'].fillna(-1)
    df['hg19chrom'] = df['hg19chrom'].fillna(-1)
    df['hg18start'] = df['hg18start'].fillna(-1)
    df['hg19start'] = df['hg19start'].fillna(-1)
    df['hg18end'] = df['hg18end'].fillna(-1)
    df['hg19end'] = df['hg19end'].fillna(-1)
    df['hg18chrom'] = df['hg18chrom'].astype(int)
    df['hg19chrom'] = df['hg19chrom'].astype(int)
    df['hg18start'] = df['hg18start'].astype(int)
    df['hg19start'] = df['hg19start'].astype(int)
    df['hg18end'] = df['hg18end'].astype(int)
    df['hg19end'] = df['hg19end'].astype(int)

    conn = sqlite3.connect(genedb_path)  # open connection

    # save tsv to sqlite3 database
    psql.write_frame(df,  # pandas dataframe
                     'nucleotide',  # table name
                     con=conn,  # connection
                     flavor='sqlite',  # use sqlite
                     if_exists='replace')  # drop table if exists

    conn.close()  # close


def main():
    """Concatenates all the mutation data from tab delmited files in
    the cosmic directory. Next, saves the results to a sqlite db."""
    # get input/output configurations
    in_opts = _utils.get_input_config('input')
    cosmic_dir = in_opts['cosmic_dir']
    out_opts = _utils.get_output_config('gene_tsv')
    out_path = out_opts['gene_tsv']
    db_opts = _utils.get_db_config('genes')
    out_db = db_opts['db']

    # save info into a txt file and sqlite3 database
    concatenate_genes(out_path, cosmic_dir)  # concatenate all gene files
    save_db(out_path, out_db)  # save concatenate file to sqlite database
