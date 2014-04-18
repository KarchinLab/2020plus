import utils.python.util as _utils
import pandas as pd
import pandas.io.sql as psql
import sqlite3


def filter_hypermutators(hypermutator_count, conn, db_path=''):
    """Query database to find hypermutator samples so they can
    be excluded from further analysis.

    **Parameters**

    hypermutator_count : int
        samples with mutation counts below this number are allowed
    conn : db connection
        database connection
    db_path : str
        if using non-config defined db, specify the db path
    """
    sql = ("SELECT *"
          " FROM maf_mutation"
          " WHERE Tumor_Sample in ("
          "     SELECT y.Tumor_Sample"
          "     FROM ("
          "         SELECT x.Tumor_Sample, SUM(x.mut_indicator) as MutationCounts"
          "         FROM ( "
          "             SELECT Tumor_Sample, 1 as mut_indicator"
          "             FROM maf_mutation"
          "         ) x "
          "         GROUP BY Tumor_Sample"
          "     ) y"
          "     WHERE y.MutationCounts<%d"
          " )" % hypermutator_count)

    df = psql.frame_query(sql, conn)  # get non hypermutator mutations

    _utils.drop_table('maf_mutation', db_path, kind='sqlite')

    psql.write_frame(df,
                     'maf_mutation',
                     conn,
                     flavor='sqlite',
                     if_exists='replace')


def save_db(maf_path, db_path, hypermutator_count):
    # merge all data frames together with the first
    # data frames given priority over later data frames
    df_cols = ['Gene_Symbol', 'Tumor_Sample', 'Tumor_Type', 'Chromosome',
               'Start_Position', 'End_Position', 'Variant_Classification',
               'Reference_Allele', 'Tumor_Allele', 'Protein_Change']
    df = pd.DataFrame(columns=df_cols)
    for single_maf in maf_path.split(','):
        tmp_df = pd.read_csv(single_maf, sep='\t')
        samp_names = set(df['Tumor_Sample'].tolist())
        tmp_df = tmp_df[tmp_df['Tumor_Sample'].apply(lambda x: x not in samp_names)]
        df = pd.concat([df, tmp_df])
        #df = df.groupby(df.index).first()  # use only first appearance of sample name

    _utils.drop_table('maf_mutation', db_path, kind='sqlite')
    conn = sqlite3.connect(db_path)  # open connection

    # save tsv to sqlite3 database
    psql.write_frame(df,  # pandas dataframe
                     'maf_mutation',  # table name
                     con=conn,  # connection
                     flavor='sqlite',  # use sqlite
                     if_exists='replace')  # drop table if exists

    # filter hypermutator samples
    filter_hypermutators(hypermutator_count, conn, db_path)


def create_empty_maf_mutation_table(db_path):
    # specify columns in table
    cols_of_interest = ['Gene_Symbol', 'Tumor_Sample', 'Tumor_Type',
                        'Chromosome', 'Start_Position', 'End_Position',
                        'Variant_Classification', 'Reference_Allele',
                        'Tumor_Allele', 'Protein_Change']
    data_type = ['TEXT', 'TEXT', 'TEXT', 'INTEGER',
                 'INTEGER', 'INTEGER', 'TEXT', 'TEXT',
                 'TEXT', 'TEXT']
    _utils.creat_empty_table('maf_mutation', db_path,
                             cols_of_interest, data_type)

def main(maf_path, db_path, hypermutator_count):
    # get db info
    db_opts = _utils.get_db_config('2020plus')
    out_db = db_opts['db']
    out_db = db_path if db_path else out_db

    # update databse maf_mutation table
    if maf_path:
        # add to database if they specify a MAF file
        save_db(maf_path, out_db, hypermutator_count)
    else:
        # else just create an empty maf_mutation table
        create_empty_maf_mutation_table(out_db)
