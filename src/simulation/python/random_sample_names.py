import numpy as np
import pandas as pd
import pandas.io.sql as psql
import src.data_analysis.python.feature_matrix as fmat
import src.data_analysis.python.position_entropy as pentropy
import src.features.python.features as feat
from src.utils.python.amino_acid import AminoAcid
import logging

logger = logging.getLogger(__name__)

class RandomSampleNames(object):

    def __init__(self, sub_sample,
                 num_iter,
                 db_conn,
                 table_name='mutation',
                 col_name='Tumor_Sample'):
        self.db_conn = db_conn
        self.set_sub_sample(sub_sample)
        self.set_num_iter(num_iter)
        self.TABLE_NAME = table_name
        self.COLUMN_NAME = col_name
        self.set_df(None)

    def dataframe_generator(self):
        """Generate subsampled data frames according to the sub_sample
        and num_iter arguments.
        """
        #n = int(self.sub_sample * self.total_count)  # number of counts to sample

        for i in range(self.num_iter):
            logger.info('Feature generation: Sub-sample rate={0}, Iteration={1} . . .'.format(self.sub_sample, i))
            # get sample names to be used
            prng = np.random.RandomState()
            prng.shuffle(self.sample_names)

            samps_of_interest = set(self.sample_names[:int(self.num_sample_names*self.sub_sample)])

            # get data from those sample names
            samp_flag = self.df[self.COLUMN_NAME].apply(lambda x: x in samps_of_interest)
            ixs = samp_flag[samp_flag==True].index
            self.current_df = self.df.ix[ixs].copy()

            # process features
            feat_list = fmat.generate_feature_matrix(self.current_df, 2)
            headers = feat_list.pop(0)  # remove header row
            feat_df = pd.DataFrame(feat_list, columns=headers)  # convert to data frame
            proc_feat_df = feat.process_features(feat_df, 0)
            miss_ent_df = pentropy.missense_position_entropy(self.current_df[['Gene', 'AminoAcid']])
            mut_ent_df = pentropy.mutation_position_entropy(self.current_df[['Gene', 'AminoAcid']])

            # encorporate entropy features
            proc_feat_df['mutation position entropy'] = mut_ent_df['mutation position entropy']
            proc_feat_df['pct of uniform mutation entropy'] = mut_ent_df['pct of uniform mutation entropy']
            proc_feat_df['missense position entropy'] = miss_ent_df['missense position entropy']
            proc_feat_df['pct of uniform missense entropy'] = miss_ent_df['pct of uniform missense entropy']

            # drop gene name column
            # proc_feat_df = proc_feat_df.drop('gene', axis=1)

            logger.info('Finished feature generation: Sub-sample rate={0}, Iteration={1}'.format(self.sub_sample, i))
            yield proc_feat_df

    def set_sub_sample(self, sub_sample):
        """Set the fraction of the original total mutations to actually sample.

        Sampling is done without replacement.

        **Parameters**

        sub_sample : float
            0 < sub_sample <= 1.0
        """
        if 0 <= sub_sample <= 1:
            self.sub_sample = sub_sample
        else:
            raise ValueError('Subsample should be between zero and one.')

    def set_num_iter(self, num_iter):
        """Set number of times to sample w/o replacement.

        **Parameters**

        num_iter : int
            do sample w/o replacement, num_iter number of times
        """
        if iter > 0:
            self.num_iter = num_iter
        else:
            raise ValueError('Number of iterations should be positive.')

    def set_df(self, df):
        if df:
            self.df = df
        else:
            sql = ("SELECT Gene, Protein_Change as AminoAcid, "
                   "       DNA_Change as Nucleotide, "
                   "       Variant_Classification, "
                   "       Tumor_Sample "
                   "FROM {0}".format(self.TABLE_NAME))
            self.df = psql.frame_query(sql, con=self.db_conn)
        self.df['is_non_silent'] = self.df['AminoAcid'].apply(lambda x: int(AminoAcid(x).is_non_silent))
        self.sample_names = self.df[self.COLUMN_NAME].unique()
        self.num_sample_names = len(self.sample_names)
        self.total_count = len(self.df)
