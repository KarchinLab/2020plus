import src.utils.python.util as _utils
import src.features.python.feature_utils as futils
import sqlite3
import pandas.io.sql as psql
import logging

logger = logging.getLogger(__name__)



#######################################
# Function used to process COSMIC data
#######################################

def main_cosmic(options):
    """Main function used to process COSMIC data."""
    # get configs
    in_opts = _utils.get_input_config('classifier')
    out_opts = _utils.get_output_config('features')
    count_opts = _utils.get_output_config('feature_matrix')
    # result_opts = _utils.get_input_config('result')
    db_cfg = _utils.get_db_config('2020plus')

    # get mutations
    conn = sqlite3.connect(db_cfg['db'])
    sql = ("SELECT Gene, Protein_Change as AminoAcid, "
            "       DNA_Change as Nucleotide, "
            "       Variant_Classification, "
            "       Tumor_Sample, Tumor_Type "
            "FROM mutations")
    mut_df = psql.frame_query(sql, con=conn)
    conn.close()

    # get features for classification
    all_features = futils.generate_features(mut_df, options)

    # save features to text file
    cols = all_features.columns.tolist()
    new_order = ['gene'] + cols[:cols.index('gene')] + cols[cols.index('gene')+1:]
    all_features = all_features[new_order]  # make the gene name the first column
    out_path = _utils.save_dir + in_opts['gene_features'] if not options['output'] else options['output']
    all_features.to_csv(out_path, sep='\t', index=False)
