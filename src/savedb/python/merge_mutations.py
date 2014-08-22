import src.utils.python.util as _utils
import sqlite3


def main(db_path):
    db_opts = _utils.get_db_config('2020plus')
    out_db = db_opts['db']
    out_db = db_path if db_path else out_db

    # drop table if exists
    # _utils.drop_table('mutation', out_db, kind='sqlite')

    cols_of_interest = ['Gene', 'Tumor_Sample', 'Tumor_Type',
                        'Chromosome', 'Start_Position',
                        'End_Position', 'Variant_Classification',
                        'Reference_Allele', 'Tumor_Allele',
                        'Protein_Change', 'DNA_Change']
    data_type = ['TEXT', 'TEXT', 'TEXT', 'TEXT', 'INT',
                 'INT', 'TEXT', 'TEXT', 'TEXT', 'TEXT', 'TEXT']
    _utils.create_empty_table('mutation', out_db,
                              cols_of_interest, data_type)

    conn = sqlite3.connect(out_db)  # open connection
    cur = conn.cursor()
    #col_info_list = [' '.join(x) for x in zip(cols_of_interest, data_type)]
    #col_info_str = ', '.join(col_info_list)
    #sql = "CREATE TABLE mutation({0});".format(col_info_str)
    #cur.execute(sql)
    #conn.commit()

    maf_mutation_cols = ['Gene_Symbol'] + cols_of_interest[1:-1]
    sql = ("INSERT INTO mutation ({0}) "
           "    SELECT {1}"
           "    FROM maf_mutation".format(', '.join(cols_of_interest[:-1]),
                                          ', '.join(maf_mutation_cols)))
    cur.execute(sql)
    conn.commit()

    cols_of_interest = cols_of_interest[:-5] + cols_of_interest[-2:] + ['Variant_Classification']
    cosmic_col_list = ['Gene', 'SampleName', 'PrimaryTissue', 'hg19chrom',
                       'hg19start', 'hg19end', 'AminoAcid', 'Nucleotide',
                       'Variant_Classification']
    mut_cols = ', '.join(cols_of_interest)
    cosmic_cols = ', '.join(cosmic_col_list)
    sql = ("INSERT INTO mutation({0}) "
           "    SELECT {1} "
           "    FROM cosmic_mutation cm"
           "    WHERE cm.SampleName NOT IN ("
           "        SELECT DISTINCT(m.Tumor_Sample) "
           "        FROM mutation m "
           "    ) ".format(mut_cols, cosmic_cols))
    cur.execute(sql)
    conn.commit()

    sql = ("UPDATE mutation SET DNA_Change='c.?'"
           "WHERE DNA_CHANGE IS NULL")
    cur.execute(sql)
    conn.commit()
    conn.close()
