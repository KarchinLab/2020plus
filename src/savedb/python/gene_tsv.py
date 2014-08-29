"""The gene_tsv module aggregates all of the tab delimited files
describing gene mutations that form the COSMIC_nuc database. It
saves the resulting tab delimited file into the data directory.
It subsequently makes a database from the file.

Genes are not used if:
 * mutations for alternative isoforms "_ENST..."

Rows are filtered if:
 * either the amino acid or nucleotide hgvs column is empty
 * both columns indicate unkown effect c.? and p.?
 * SomaticStatus has word 'unkown'
 * row is a part of a sample defined as a hypermutator (>500)

NOTE: This module will likely break if minor details in the
database, etc. change.
"""

import src.utils.python.util as _utils
from collections import namedtuple
import pandas as pd
import pandas.io.sql as psql
import sqlite3
import string
import os
import re


def handle_cosmic_mutation_export(df, only_genome_wide, use_unknown_status):
    # rename columns
    rename_dict = {'Gene name': 'Gene',
                   'Genome-wide screen': 'GenomeWideScreen',
                   'Sample name': 'SampleName',
                   'ID_sample': 'COSMICSampleID',
                   'Mutation AA': 'AminoAcid',
                   'Mutation CDS': 'Nucleotide',
                   'Primary site': 'PrimaryTissue',
                   'Mutation somatic status': 'SomaticStatus',
                   'Mutation GRCh37 strand': 'strand',
                   'Mutation GRCh37 genome position': 'genome_pos',
                   'Mutation Description': 'Variant_Classification'}
    df = df.rename(columns=rename_dict)

    # parse genome coordinates
    is_pos = df['genome_pos'].notnull()
    genome_pos = df[is_pos]['genome_pos']
    mychr = genome_pos.apply(lambda x: x.split(':')[0])
    mystart = genome_pos.apply(lambda x: x.split(':')[1].split('-')[0])
    myend = genome_pos.apply(lambda x: x.split(':')[1].split('-')[1])

    # assign columns to data frame
    df['hg19chrom'] = ''
    df[is_pos]['hg19chrom'] = mychr
    df['hg19start'] = -1
    df[is_pos]['hg19start'] = mystart
    df['hg19end'] = -1
    df[is_pos]['hg19end'] = myend

    # if only genome wide screen flag is true then remove all non
    # genome wide data
    if only_genome_wide:
        df = df[df['GenomeWideScreen']=='y']

    # filter out genes with _ENST in name
    df = df[df['Gene'].apply(lambda x: '_ENST' not in x)]

    # filter out designated germline variants
    df = df[df['SomaticStatus'].apply(lambda x: 'germline' not in x.lower())]

    if not use_unknown_status:
        # filter out unknown somatic status mutations
        df = df[df['SomaticStatus'].apply(lambda x: 'unknown' not in x.lower())]

    # correct mutation name
    mut_type_dict = {'Substitution - Missense': 'Missense_Mutation',
                     'Substitution - coding silent': 'Silent',
                     'Substitution - Nonsense': 'Nonsense_Mutation',
                     'Frameshift': 'Frame_Shift_Indel',
                     'Deletion - Frameshift': 'Frame_Shift_Indel',
                     'Insertion - Frameshift': 'Frame_Shift_Indel',
                     'Complex - frameshift': 'Frame_Shift_Indel',
                     'Deletion - In frame': 'In_Frame_Indel',
                     'Insertion - In frame': 'In_Frame_Indel',
                     'Complex - insertion inframe': 'In_Frame_Indel',
                     'Complex - deletion inframe': 'In_Frame_Indel'}
    df['Variant_Classification'] = df['Variant_Classification'].apply(lambda x: mut_type_dict[x] if (x in mut_type_dict) else x)

    col_order = ['Gene', 'SampleName', 'COSMICSampleID', 'AminoAcid',
                 'Nucleotide', 'PrimaryTissue', 'SomaticStatus',
                 'strand', 'hg19chrom', 'hg19start', 'hg19end',
                 'Variant_Classification']
    return df[col_order]


def parse_sample_name(sample):
    if sample.startswith('TCGA'):
        sample = re.sub('-[A-Za-z0-9]+$', '', sample)
        return sample
    else:
        return sample


def skip_header(file_handle, skip_rows=8):
    """Skips the first "skip_row" lines of the file.

    To skip rows the next() method is called repeatedly.
    This skips the non-table part of the tab delimited gene files.

    Parameters
    ----------
    file_handle : file
        file object of tab delim file
    skip_rows : int, Default=8
        number of lines to skip
    """
    for i in range(skip_rows):
        file_handle.next()
    return file_handle


def concatenate_genes(out_path, cosmic_dir):
    """Saves all information in tab delimited files from genes into
    a single file in the data directory.

    Exceptions
    * Skips "genes" that look like Ensemble Transcripts ("_ENST")
    * Skips rows with "unkown" in somatic status
    * Skips rows that are blank at either AA or nucleotide mutations
    * Skips rows where mutations are both "?" for AA and DNA

    Parameters
    ----------
    out_path : str
        path to save concatenated file
    cosmic_dir : str
        base directory for gene tsv files
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
                            LineTuple = namedtuple('LineTuple', header[1:])  # create namedtuple type
                            FIRST_FLAG = False

                        gene_name = file_name[:-4]  # file name before ".tsv"
                        handle = skip_header(handle)  # skip beginning lines
                        # add data
                        for line in handle:
                            #split_line = line.split('\t')
                            split_tuple = LineTuple(*line.split('\t'))  # split
                            if split_tuple.AminoAcid and split_tuple.Nucleotide:
                                # if line designates a mutation
                                if split_tuple.AminoAcid != 'p.?' or split_tuple.Nucleotide != 'c.?':
                                    # not unknown effect for AA and amino acid
                                    if 'unknown' not in split_tuple.SomaticStatus.lower():
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
                        #split_line = line.split('\t')
                        split_tuple = LineTuple(*line.split('\t'))  # split line
                        if split_tuple.AminoAcid and split_tuple.Nucleotide:
                            # if line designates a mutation
                            if split_tuple.AminoAcid != 'p.?' or split_tuple.Nucleotide != 'c.?':
                                # not uknown effect for both AA and nucleotide
                                if 'unknown' not in split_tuple.SomaticStatus.lower():
                                    # do not include unknown somatic status mutations
                                    mywriter.write(gene_name + "\t" + line)  # write with gene name


def filter_hypermutators(hypermutator_count, conn, db_path=''):
    """Query database to find hypermutator samples so they can
    be excluded from further analysis.

    Parameters
    ----------
    hypermutator_count : int
        samples with mutation counts below this number are allowed
    conn : db connection
        database connection
    db_path : str
        if using non-config defined db, specify the db path
    """
    sql = ("SELECT *"
          " FROM cosmic_mutation"
          " WHERE COSMICSampleID in ("
          "     SELECT y.COSMICSampleID"
          "     FROM ("
          "         SELECT x.COSMICSampleID, SUM(x.mut_indicator) as MutationCounts"
          "         FROM ( "
          "             SELECT COSMICSampleID, 1 as mut_indicator"
          "             FROM cosmic_mutation"
          "         ) x "
          "         GROUP BY COSMICSampleID"
          "     ) y"
          "     WHERE y.MutationCounts<%d"
          " )" % hypermutator_count)

    df = psql.frame_query(sql, conn)  # get non hypermutator mutations

    _utils.drop_table('cosmic_mutation', db_path, kind='sqlite')

    psql.write_frame(df,
                     'cosmic_mutation',
                     conn,
                     flavor='sqlite',
                     if_exists='replace')


def save_db(hypermutator_ct,
            gene_tsv_path,
            genedb_path,
            is_genes_tgz=False,
            only_genome_wide=True,
            use_unknown_status=False):
    """Saves tab delim gene mutation file to a sqlite3 db.

    NOTE: Uses pandas to store all contents in memory and then
    saves to sqlite db. This may cause large memory usage.

    Parameters
    ----------
    hypermutator_ct : int
        filter for overly mutated samples
    gene_tsv_path : str
        path to tab delim file containing all gene mutations
    cnv_tsv_path : str
        path to tab delim file containing cosmic cnv mutations
    cell_line_path : str
        path to cosmic cell line project file
    genedb_pah : str
        path to sqlite3 db
    """
    # read data
    df = pd.read_csv(gene_tsv_path, sep='\t')
    # cnv_df = pd.read_csv(cnv_tsv_path, sep=r'\t|:|\.\.')

    # filter out cell line samples
    #if cell_line_path:
        #cell_line_df = pd.read_csv(cell_line_path, sep='\t')
        #cell_line_sample_names = set(cell_line_df['Sample name'].tolist())
    #else:
        #cell_line_sample_names = set([])

    # skip this if COSMIC not used
    #df = df[df['SampleName'].apply(lambda x: x not in cell_line_sample_names)]

    if is_genes_tgz:
        # fix sample names so they match with external data
        df['SampleName'] = df['SampleName'].apply(parse_sample_name)

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
    else:
        df = handle_cosmic_mutation_export(df, only_genome_wide, use_unknown_status)

    # drop table if already exists
    _utils.drop_table('cosmic_mutation', genedb_path, kind='sqlite')
    # _utils.drop_table('cosmic_cnv', genedb_path, kind='sqlite')

    conn = sqlite3.connect(genedb_path)  # open connection

    # save tsv to sqlite3 database
    psql.write_frame(df,  # pandas dataframe
                     'cosmic_mutation',  # table name
                     con=conn,  # connection
                     flavor='sqlite',  # use sqlite
                     if_exists='replace')  # drop table if exists
    #psql.write_frame(cnv_df,  # pandas dataframe
                     #'cosmic_cnv',  # table name
                     #con=conn,  # connection
                     #flavor='sqlite',  # use sqlite
                     #if_exists='replace')  # drop table if exists

    # drop table and re-insert data without hypermutators
    filter_hypermutators(hypermutator_ct, conn, genedb_path)

    conn.close()  # close


def create_empty_cosmic_mutation_table(db_path):
    # construct column/column types
    mycols = ['Gene', 'SampleName', 'COSMICSampleID', 'AminoAcid',
              'Nucleotide', 'PrimaryTissue', 'Tissuesubtype1',
              'Tissuesubtype2', 'Histology', 'Histologysubtype1',
              'Histologysubtype2', 'PubmedID', 'studies',
              'MutationID', 'SomaticStatus', 'SampleSource',
              'Zygosity', 'hg18chrom', 'hg18start', 'hg18end',
              'hg19chrom', 'hg19start', 'hg19end', 'Variant_Classification']
    coltypes = ['TEXT']*len(mycols)
    coltypes[-1] = 'INTEGER'
    coltypes[-2] = 'INTEGER'
    coltypes[-4] = 'INTEGER'
    coltypes[-5] = 'INTEGER'

    # create empty cosmic_mutation table
    _utils.create_empty_table('cosmic_mutation', db_path,
                              mycols, coltypes)


def main(hypermutator_count, mut_path, db_path,
         no_cosmic_flag, opts):
    """Concatenates all the mutation data from tab delmited files in
    the cosmic directory. Next, saves the results to a sqlite db.

    Parameters
    ----------
    hypermutator_count : int
        remove samples with too many mutations
    mut_path : str
        Either path to directory containing contents of COSMIC's
        genes.tgz file or decompressed CosmicMutantExport.tsv.
        If empty string, just use path from config file.
    db_path : str
        path to save sqlite database. If string is empty,
        use path from config.
    no_cosmic_flag : bool
        indicates not to use cosmic mutations
    """
    # get input/output configurations
    in_opts = _utils.get_input_config('input')
    cosmic_path = in_opts['cosmic_path']
    out_opts = _utils.get_output_config('gene_tsv')
    out_path = out_opts['gene_tsv']
    cnv_path = out_opts['cnv_tsv']
    db_opts = _utils.get_db_config('2020plus')
    out_db = db_opts['db']

    # check if user specifies non standard db path
    out_db = db_path if db_path else out_db

    # save info into a txt file and sqlite3 database
    if not no_cosmic_flag:
        cosmic_path = mut_path if mut_path else cosmic_path
        if os.path.isdir(cosmic_path):
            # concatenate all gene files
            concatenate_genes(out_path, cosmic_path)
            # save database
            save_db(hypermutator_count, out_path, out_db,
                    is_genes_tgz=True,
                    only_genome_wide=opts['only_genome_wide'],
                    use_unknown_status=opts['use_unknown_status'])
        elif os.path.isfile(cosmic_path):
            # save database
            save_db(hypermutator_count, cosmic_path, out_db,
                    is_genes_tgz=False,
                    only_genome_wide=opts['only_genome_wide'],
                    use_unknown_status=opts['use_unknown_status'])
        else:
            raise ValueError('Please specify a vlid path to COSMIC data')
    else:
        # create an empty table if cosmic not wanted
        create_empty_cosmic_mutation_table(out_db)
