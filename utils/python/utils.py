import pandas as pd
import logging


def read_aa_properties(file_path):
    """Read aa property counts from the data_analysis/results folder.

    Args:
        file_path (str): path to aa_change.properties.txt

    Returns:
        DataFrame. contains mutation counts for amino acid chemical properties
    """
    logger = logging.getLogger(name=__name__)
    logger.info('reading in %s ...' % file_path)
    df = pd.read_csv(file_path, sep='\t')  # read file
    df = df.set_index('initial_prop')  # set rows as initial property
    logger.info('finished reading file.')
    return df
