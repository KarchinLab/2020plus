"""The gene_tsv module aggregates all of the tab delimited files
describing gene mutations that form the COSMIC_nuc database. It
saves the resulting tab delimited file into the data directory.
It subsequently makes a database from the file.

NOTES
 * only lines with either a nucleotide or amino acid mutation
   are saved.
 * this module in basically just a script
 * it will likely break if paths/conventtions change
 * it is used because the gene name needs to be added to the
   concatenated file. Otherwise the "cat" command could have
   been used.
"""

import utils.python.util as _utils
import pandas as pd
import pandas.io.sql as psql
import numpy as np
import sqlite3
import string
import os


def skip_header_next(file_handle, skip_rows=8):
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


def skip_header_readline(file_handle, skip_rows=8):
    """Skips the first "skip_row" lines of the file.

    To skip rows the readline() method is called repeatedly.
    This skips the non-table part of the tab delimited gene files.

    Args
        file_handle (file): file object of tab delim file

    Kwargs
        skip_rows (int): number of lines to skip
    """
    for i in range(skip_rows):
        file_handle.readline()
    return file_handle


def concatenate_genes(out_path, cosmic_dir):
    """Saves all information in tab delimited files from genes into
    a single file in the data directory.
    """
    # concatenate gene tab delim files
    with open(out_path, 'wb') as mywriter:
        # iterate through 'A'...'Z' directories
        FIRST_FLAG = True  # flag for including header
        for letter in string.ascii_uppercase:
            for file_name in os.listdir(cosmic_dir + letter):
                if file_name.endswith('.tsv') and '_ENST' not in file_name:
                    # only use tab delim files that are
                    # not alternative isoforms
                    with open(cosmic_dir + letter + "/" + file_name) as handle:
                        gene_name = file_name[:-4]  # file name before ".tsv"
                        if FIRST_FLAG:
                            # first file
                            handle = skip_header_readline(handle, skip_rows=7)
                            line = handle.readline()
                            mywriter.write('genes\t' + line)
                            FIRST_FLAG = False  # only include header once
                            line = handle.readline()
                            while line:
                                split_line = line.split('\t')
                                if split_line[2] or split_line[3]:
                                    # if line designates a mutation
                                    mywriter.write(gene_name + "\t" + line)
                                line = handle.readline()
                        else:
                            handle = skip_header_next(handle)  # skip beginning lines
                            # add data
                            for line in handle:
                                split_line = line.split('\t')
                                if split_line[2] or split_line[3]:
                                    # if line designates a mutation
                                    mywriter.write(gene_name + "\t" + line)  # write with gene name


        # iterate through genes in the numeric directory
        numeric_dir = cosmic_dir + '0-9/'
        for file_name in os.listdir(numeric_dir):
            if file_name.endswith('.tsv') and '_ENST' not in file_name:
                # only use tab delim files that are
                # not alternative isoforms
                with open(numeric_dir + file_name) as handle:
                    gene_name = file_name[:-4]  # file name before ".tsv"
                    handle = skip_header_next(handle)  # skip beginning lines
                    for line in handle:
                        split_line = line.split('\t')
                        if split_line[2] or split_line[3]:
                            # if line designates a mutation
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
    df['Chromosome NCBI36'] = df['Chromosome NCBI36'].fillna(-1)
    df['Chromosome GRCh37'] = df['Chromosome GRCh37'].fillna(-1)
    df['Genome Start NCBI36'] = df['Genome Start NCBI36'].fillna(-1)
    df['Genome Start GRCh37'] = df['Genome Start GRCh37'].fillna(-1)
    df['Genome Stop NCBI36'] = df['Genome Stop NCBI36'].fillna(-1)
    df['Genome Stop GRCh37'] = df['Genome Stop GRCh37'].fillna(-1)
    df['Chromosome NCBI36'] = df['Chromosome NCBI36'].astype(int)
    df['Chromosome GRCh37'] = df['Chromosome GRCh37'].astype(int)
    df['Genome Start NCBI36'] = df['Genome Start NCBI36'].astype(int)
    df['Genome Start GRCh37'] = df['Genome Start GRCh37'].astype(int)
    df['Genome Stop NCBI36'] = df['Genome Stop NCBI36'].astype(int)
    df['Genome Stop GRCh37'] = df['Genome Stop GRCh37'].astype(int)

    conn = sqlite3.connect(genedb_path)  # open connection

    # save tsv to sqlite3 database
    psql.write_frame(df,  # pandas dataframe
                     'mutations',  # table name
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
    out_db = out_opts['gene_db']

    # save info into a txt file and sqlite3 database
    concatenate_genes(out_path, cosmic_dir)  # concatenate all gene files
    save_db(out_path, out_db)  # save concatenate file to sqlite database
