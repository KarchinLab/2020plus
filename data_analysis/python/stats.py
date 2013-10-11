from db_utils.cosmic_db import get_cosmic_db

def count_mutations(cursor):
    """Count the number of entries"""
    cursor.execute("""SELECT COUNT(COSMICSampleID)
                   FROM `nucleotide`""")
    return cursor.fetchone()[0]  # COUNT query returns a tuple

def main():
    conn = get_cosmic_db()
    cursor = conn.cursor()
    mycount = count_mutations(cursor)  # count mutations
    print mycount
    conn.close()

if __name__=="__main__":
    main()
