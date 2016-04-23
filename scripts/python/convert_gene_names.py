import pandas as pd
import numpy as np
import csv
import argparse
import IPython


already_converted = {}


def parse_arguments():
    info = 'Convert gene names to approved HUGO symbol'
    parser = argparse.ArgumentParser(description=info)
    help_str = 'Path to HUGO name file'
    parser.add_argument('-hugo', '--hugo',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Path to input file'
    parser.add_argument('-i', '--input',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Column header name or column number to be converted'
    parser.add_argument('-c', '--column',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Prevent duplication of genes when converting'
    parser.add_argument('-no-duplication', '--no-duplication',
                        default=False, action='store_true',
                        help=help_str)
    help_str = 'New output with converted gene names'
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def convert_genes(gname, approved_symbols, symbol_dict, synonym_dict):
    if gname in approved_symbols:
        return gname
    elif gname in symbol_dict:
        symb = symbol_dict[gname]
        already_converted.setdefault(symb, gname)
        return symb
    elif gname in synonym_dict:
        syn = synonym_dict[gname]
        already_converted.setdefault(syn, gname)
        return syn
    else:
        return gname


def main(opts):
    # make gene conversion dictionaries
    hugo_df = pd.read_csv(opts['hugo'], sep='\t')
    hugo_approved = hugo_df['Approved Symbol']
    prev_symbols = hugo_df['Previous Symbols'].str.split(', ')
    synonyms = hugo_df['Synonyms'].str.split(', ')
    symbols2hugo = {}
    synonyms2hugo = {}
    for ix in range(len(hugo_df)):
        # add previous symbol to conversion dictionary
        if prev_symbols[ix] is not np.nan:
            for psymbol in prev_symbols[ix]:
                symbols2hugo[psymbol] = hugo_approved[ix]

        # add synonyms to conversion dictionary
        if synonyms[ix] is not np.nan:
            for syn in synonyms[ix]:
                synonyms2hugo[syn] = hugo_approved[ix]
    # fix stupid KMT2D/KMT2B MLL2/MLL4 problem
    symbols2hugo['MLL4'] = 'KMT2B'
    synonyms2hugo['MLL4'] = 'KMT2B'

    # figure out which column to convert
    if opts['column'].isdigit():
        col_num = int(opts['column'])
        output = []

        # read in input
        with open(opts['input']) as handle:
            approved_hugos = set(hugo_approved.tolist())
            for line in handle:
                line_split = line.strip('\n').split('\t')
                tmp = convert_genes(line_split[col_num],
                                    approved_hugos,
                                    symbols2hugo, synonyms2hugo)
                # add to output
                if tmp == line_split[col_num]:
                    line_split[col_num] = tmp
                    output.append(line_split)
                elif line_split[col_num] == already_converted[tmp]:
                    line_split[col_num] = tmp
                    output.append(line_split)

        # write output
        with open(opts['output'], 'w') as handle:
            csv.writer(handle, delimiter='\t', lineterminator='\n').writerows(output)
    else:
        input_df = pd.read_csv(opts['input'], sep='\t')
        myargs = (set(hugo_approved.tolist()), symbols2hugo, synonyms2hugo)
        input_df[opts['column']] = input_df[opts['column']].apply(convert_genes, args=myargs)
        input_df.to_csv(opts['output'], sep='\t', index=False)


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
