import numpy as np
import features.python.features as features

class Bootstrap(object):

    def __init__(self, df, subsample=1.0, num_iter=10):
        # initialize variables
        self.set_df(df)
        self.set_subsample(subsample)
        self.set_num_iter(num_iter)

        # set up data for multinomial sampling
        self.nonzero_indices = np.nonzero(self.df.values)
        counts = self.df.values[self.nonzero_indices]
        self.total_count = np.sum(counts)
        self.prob = counts.astype(float) / self.total_count

    def dataframe_generator(self):
        """Generate subsampled data frames according to the subsample
        and num_samples arguments.
        """
        tmp_df = self.df.copy()  # copy df so it is not mutable
        n = int(self.subsample * self.total_count)  # number of counts to sample
        for i in range(self.num_iter):
            prng = np.random.RandomState()
            multinomial_sample = prng.multinomial(n,  # total counts for multinomial
                                                  self.prob)  # probability
            tmp_df.values[self.nonzero_indices] = multinomial_sample
            tmp_df = features.process_features(tmp_df, 0)
            yield tmp_df

    def set_subsample(self, subsample):
        """Set the fraction of the original total count to actually sample.

        Sampling is done with replacement. If subsample == 1.0 then it is a
        traditional bootstrap procedure. However, subsamples < 1.0 will
        identify dependencies on database size.

        Args:
          | subsample (float): 0 < subsample <= 1.0
        """
        if subsample > 0:
            self.subsample = subsample
        else:
            raise ValueError('Subsample should be positive.')

    def set_num_iter(self, num_iter):
        """Set number of times to sample from multinomial distribution.

        Args:
          | num_samples (int): get multinomial sample, num_samples number of times
        """
        if num_iter > 0:
            self.num_iter = num_iter
        else:
            raise ValueError('Number of iterations should be positive.')

    def set_df(self, df):
        """Set data frame which contains the counts for mutational types."""
        if 'gene' in df.columns:
            df = df.set_index('gene')
        self.df = df
