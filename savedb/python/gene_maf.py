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
    df = pd.read_csv(maf_path, sep='\t')  # read data
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

def main(maf_path, db_path, hypermutator_count):
    db_opts = _utils.get_db_config('champ')
    out_db = db_opts['db']
    out_db = db_path if db_path else out_db
    save_db(maf_path, out_db, hypermutator_count)
