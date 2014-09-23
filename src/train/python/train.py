from src.classify.python.r_random_forest_clf import RRandomForest
import src.utils.python.util as _utils
import pandas as pd
import logging

logger = logging.getLogger(__name__)

def main(cli_opts):
    cfg_opts = _utils.get_output_config('classifier')
    in_opts = _utils.get_input_config('classifier')
    minimum_ct = cli_opts['min_count']

    # read in features
    df = pd.read_csv(_utils.save_dir + in_opts['gene_feature'],
                     sep='\t', index_col=0)

    logger.info('Training R\'s Random forest . . .')
    rrclf = RRandomForest(df,
                          other_sample_ratio=cli_opts['other_ratio'],
                          driver_sample=cli_opts['driver_rate'],
                          ntrees=cli_opts['ntrees'],
                          min_ct=minimum_ct)
    rrclf.train()  # train on entire data
    logger.info('Finished training.')
    logger.info('Saving classifier to . . .')
    rrclf.clf.save(cli_opts['output'])
    logger.info('Finished saving classifier.')


if __name__ == "__main__":
    main()
