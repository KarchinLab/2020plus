import MySQLdb

def get_cosmic_db():
    """Return a connection to the cosmic database."""
    conn = MySQLdb.connect(host='karchin-db01.icm.jhu.edu',
                           user='collin',
                           passwd='Bgt5r4r5',
                           db='COSMIC_nuc')
    return conn
