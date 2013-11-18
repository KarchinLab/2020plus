from vogelstein_classifier import VogelsteinClassifier
import utils.python.util as _utils
import pandas as pd


def main():
    cfg_opts = _utils.get_output_config('classifier')
    input_opts = _utils.get_output_config('stats')

    vclf = VogelsteinClassifier()  # 20/20 rule classifier
    df = pd.read_csv(_utils.result_dir + input_opts['gene_design_matrix'],
                     sep='\t', index_col=0)
    df['total'] = df.T.sum()
    df['curated class'] = [_utils.classify_gene(gene)
                           for gene in df.index.tolist()]
    input_list = ((row['recurrent missense'],
                   row['frame shift'] + row['nonsense'],
                   row['total'])
                  for i, row in df.iterrows())
    pred_classes = vclf.predict_list(input_list)
    df['2020 class'] = pred_classes
    df.to_csv(_utils.clf_result_dir + cfg_opts['vogelstein_predictions'],
              sep='\t')


if __name__ == "__main__":
    main()
