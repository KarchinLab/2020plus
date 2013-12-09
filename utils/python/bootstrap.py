import numpy as np

class Bootstrap(object):

    def __init__(self, df, subsample=1.0, num_samples=10):
        # initialize variables
        self.set_df(df)
        self.set_subsample(subsample)
        self.set_num_samples(num_samples)

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
        multinomial_sample = np.random.multinomial(n,  # total counts for multinomial
                                                   self.prob,  # probability
                                                   self.num_samples)  # number of samples
        for i in range(self.n_iter):
            tmp_df.values[self.nonzero_indices] = multinomial_sample[i]
            yield tmp_df

    def set_subsample(self, subsample):
        """Set the fraction of the original total count to actually sample.

        Sampling is done with replacement. If subsample == 1.0 then it is a
        traditional bootstrap procedure. However, subsamples < 1.0 will
        identify dependencies on database size.

        Args:
          | subsample (float): 0 < subsample <= 1.0
        """
        if 0 < subsample <= 1:
            self.subsample = subsample
        else:
            raise ValueError('Subsample should be between zero and one.')

    def set_num_samples(self, num_samples):
        """Set number of times to sample from multinomial distribution.

        Args:
          | num_samples (int): get multinomial sample, num_samples number of times
        """
        if num_samples > 0:
            self.num_samples = num_samples
        else:
            raise ValueError('Number of iterations should be positive.')

    def set_df(self, df):
        """Set data frame which contains the counts for mutational types."""
        self.df = df

