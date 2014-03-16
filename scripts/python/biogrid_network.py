import networkx as nx
import pandas as pd
import argparse
import csv


def parse_arguments():
    description = 'Creates a biogrid interaction network.'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-b', '--biogrid',
                        type=str, required=True,
                        help='Biogrid text file in tab2 format')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output path to save features')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    df = pd.read_csv(opts['biogrid'], sep='\t')
    interact_df = df[['Official Symbol Interactor A',
                      'Official Symbol Interactor B']]
    interact_genes = interact_df.dropna().values.tolist()
    G = nx.Graph()
    G.add_edges_from(map(tuple, interact_genes))
    gene_betweeness = nx.betweenness_centrality(G)
    gene_degree = G.degree()
    result = [[key, gene_betweeness[key], gene_degree[key]]
              for key in gene_degree]
    result = [['gene', 'gene_betweeness', 'gene_degree']] + result
    with open(opts['output'], 'wb') as handle:
        csv.writer(handle, delimiter='\t').writerows(result)
    # gene_communicability = nx.communicability_centrality(G)


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
