import pandas as pd
import numpy as np
import re
import argparse


def fix_hgvs(hgvs_string):
    if hgvs_string is np.nan:
        # handle na values
        return np.nan
    if 'in_frame_ins' in hgvs_string.lower():
        hgvs_string = re.sub('in_frame_ins', 'ins', hgvs_string,
                             flags=re.IGNORECASE)
    elif 'in_frame_del' in hgvs_string.lower():
        hgvs_string = re.sub('in_frame_del', 'del', hgvs_string,
                             flags=re.IGNORECASE)
    return hgvs_string


def parse_tumor_sample(tsample):
    if tsample.startswith('TCGA'):
        tsample = re.sub('-[A-Za-z0-9]+$', '', tsample)
        return tsample
    else:
        return tsample


def parse_arguments():
    parser = argparse.ArgumentParser()
    help_str = 'Data set from Davoli et al found at Elledge lab website'
    parser.add_argument('-t', '--txt',
                        type=str, required=True,
                        help=help_str)
    help_str = 'Path to save MAF file'
    parser.add_argument('-m', '--maf',
                        type=str, required=True,
                        help=help_str)
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in data set
    davoli_df = pd.read_csv(opts['txt'], sep='\t')

    # convert positions to have separate cols for chr, start, and end
    tmp_list = [re.split(':|-', x) for x in davoli_df['Genome.position.hg19'].tolist()]
    col_names = ['Chromosome', 'Start_Position', 'End_Position']
    tmp_df = pd.DataFrame(tmp_list, columns=col_names)
    davoli_df[col_names] = tmp_df

    # convert classification of mutation types to something
    # more compatable with the MAF file format
    type2varClass = {'Silent': 'Silent',
                     'Splice Site': 'Splice_Site',
                     'Missense': 'Missense_Mutation',
                     'Indel Frameshift': 'Frame_Shift_Indel',
                     'Nonsense': 'Nonsense_Mutation',
                     'Indel In Frame': 'In_Frame_Indel',
                     'Complex Missense': 'In_Frame_Indel',
                     'Translation_Start_Site': 'Translation_Start_Site',
                     'De_novo_Start_OutOfFrame': 'Translation_Start_Site',
                     'De_novo_Start_InFrame': 'Translation_Start_Site',
                     'Nonstop Extension': 'Nonstop_Mutation',
                     'Complex Nonsense': 'Nonsense_Mutation',
                     'Indel Nonsense': 'Nonsense_Mutation',
                     "5'UTR": "5'UTR",
                     "3'UTR": "3'UTR",
                     'Intron': 'Intron',
                     'Promoter': 'Promoter'
                     }
    davoli_df['Mutation_Type'] = davoli_df['Mutation_Type'].apply(lambda x: type2varClass[x])
    davoli_df['Tumor_Sample'] = davoli_df['Tumor_Sample'].apply(parse_tumor_sample)
    davoli_df['Protein_Change'] = davoli_df['Protein_Change'].apply(fix_hgvs)

    # fix column name of ref/new base
    davoli_df.rename(columns={'Gene': 'Gene_Symbol',
                              'Mutation_Type': 'Variant_Classification',
                              'Reference': 'Reference_Allele',
                              'Mutation': 'Tumor_Allele'},
                     inplace=True)

    cols_of_interest = ['Gene_Symbol', 'Tumor_Sample', 'Tumor_Type',
                        'Chromosome', 'Start_Position',
                        'End_Position', 'Variant_Classification',
                        'Reference_Allele', 'Tumor_Allele',
                        'Protein_Change']
    davoli_df[cols_of_interest].to_csv(opts['maf'], sep='\t', index=False)


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
