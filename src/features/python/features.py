import src.utils.python.util as _utils
import src.features.python.feature_utils as futils
import sqlite3
import pandas as pd
import pandas.io.sql as psql
import numpy as np
import os
import logging

logger = logging.getLogger(__name__)

######################################
# Function used to process prob 20/20 results
######################################

def main(opts):
    # read in config file
    in_opts = _utils.get_input_config('input')

    # read in prob 20/20 files
    count_df = pd.read_csv(opts['summary'], sep='\t')
    tsg_test_df = pd.read_csv(opts['tsg_test'], sep='\t')
    og_test_df = pd.read_csv(opts['og_test'], sep='\t')
    og_test_df = og_test_df.rename(columns={'gene':'Gene'})
    tsg_test_df = tsg_test_df.rename(columns={'gene':'Gene'})

    # make feature matrix
    feature_df = futils.process_features(count_df)
    tsg_test_cols = ['Gene', 'inactivating p-value']
    feature_df = pd.merge(feature_df, tsg_test_df[tsg_test_cols],
                          how='left', on='Gene')
    og_test_cols = ['Gene', 'entropy p-value', 'vest p-value', 'combined p-value']
    feature_df = pd.merge(feature_df, og_test_df[og_test_cols],
                          how='left', on='Gene')

    # add covariate feature columns
    if opts['covariates']:
        covar_file = opts['covariates']
    else:
        covar_file = os.path.join(_utils.proj_dir, in_opts['mutsigcv_features'])
    covar_df = pd.read_csv(covar_file, sep='\t')
    covar_cols = ['gene',
                  'expression_CCLE',
                  'replication_time',
                  'HiC_compartment',
                 ]
    covar_df = covar_df[covar_cols].rename(columns={'gene': 'Gene'})
    feature_df = pd.merge(feature_df, covar_df,
                          how='left', on='Gene')

    # add biogrid features if present
    if str(opts['biogrid']).lower() != "no":
        # set biogrid features from config path if not set by user
        if opts['biogrid']:
            biogrid_file = opts['biogrid']
        else:
            biogrid_file = os.path.join(_utils.proj_dir, in_opts['biogrid_features'])
        # read in biogrid data
        biogrid_df = pd.read_csv(biogrid_file, sep='\t')
        biogrid_df = biogrid_df.rename(columns={'gene': 'Gene'})

        # permute feature if toggled
        if opts['permute_biogrid']:
            prng = np.random.RandomState(opts['random_seed'])
            bg_feats = ['gene_degree', 'gene_betweeness']
            permute_order = prng.choice(len(biogrid_df),
                                        size=len(biogrid_df),
                                        replace=False)
            biogrid_df.loc[:,bg_feats] = biogrid_df[bg_feats].loc[permute_order].values

        # merge in biogrid features
        feature_df = pd.merge(feature_df, biogrid_df, how='left', on='Gene')
        feature_df['gene_degree'] = feature_df['gene_degree'].fillna(0)
        feature_df['gene_betweeness'] = feature_df['gene_betweeness'].fillna(0)

    # fill na values
    rename_dict = {'Gene': 'gene'}
    feature_df = feature_df.rename(columns=rename_dict)
    feature_df = feature_df.fillna(feature_df.mean())

    # setup output cols reflecting feature selection
    feature_df.to_csv(opts['output'], sep='\t', index=False)


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
