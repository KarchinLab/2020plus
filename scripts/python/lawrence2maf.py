import pandas as pd
import numpy as np
import argparse

def fix_tumor_sample(tsample):
    if 'tumor' in tsample.lower():
        return tsample.strip('-Tumor')
    else:
        return tsample


def fix_tumor_type(ttype):
    ttype_dict = {'OV': 'Ovarian',
                  'ESO': 'Esophageal Adenocarcinoma',
                  'CRC': 'Colorectal',
                  'NB': 'Neuroblastoma',
                  'GBM': 'Glioblastoma Multiforme',
                  'PRAD': 'Prostate Adenocarcinoma',
                  'UCEC': 'Endometrial Carcinoma',
                  'BLCA': 'Bladder Urothelial Carcinoma',
                  'BRCA': 'Breast Adenocarcinoma',
                  'MEL': 'Melanoma',
                  'LUSC': 'Lung Squamous Cell Carcinoma',
                  'MED': 'Medulloblastoma',
                  'HNSC': 'Head and Neck Squamous Cell Carcinoma',
                  'LUAD': 'Lung Adenocarcinoma',
                  'KIRC': 'Kidney Clear Cell Carcinoma'}
    if ttype_dict.has_key(ttype):
        return ttype_dict[ttype]
    else:
        return ttype


def fix_variant_type(var_clf):
    if 'missense' in var_clf.lower():
        return 'Missense_Mutation'
    elif 'splice_site' in var_clf.lower():
        return 'Splice_Site'
    elif 'frame_shift' in var_clf.lower():
        return 'Frame_Shift_Indel'
    elif 'in_frame' in var_clf.lower():
        return 'In_Frame_Indel'
    elif 'synonymous' in var_clf.lower():
        return 'Silent'
    elif 'nonsense' in var_clf.lower():
        return 'Nonsense_Mutation'
    elif 'silent' in var_clf.lower():
        return 'Silent'
    elif 'nonstop' in var_clf.lower():
        return 'Nonstop_Mutation'
    elif 'translation_start_site' in var_clf.lower():
        return 'Translation_Start_Site'
    return var_clf


def generate_hgvs_syntax(df):
    hgvs_list = []
    for row in df.iterrows():
        row = row[1]  # grab the pandas series object
        seq_ont = row.ix['Sequence Ontology']
        nulls_in_row = row.isnull()

        # fix bug related to pandas forcing conversion of
        # dtype to float if there exists a NaN
        if not nulls_in_row['Amino acid position']:
            pos = int(row.ix['Amino acid position'])
        else:
            pos = row.ix['Amino acid position']

        if seq_ont in ['SY', 'SL', 'SG', 'MS', 'CS']:
            ref = row.ix['Reference amino acid(s)']
            alt = row.ix['Alternate amino acid(s)']
            hgvs_string = 'p.{0}{1}{2}'.format(ref, pos, alt)
            hgvs_list.append(hgvs_string)
        elif seq_ont == 'II':
            ins_len = len(row.ix['Alternate base'][:-1]) / 3
            hgvs_string = 'p.K{0}_K{1}ins{2}'.format(pos, pos+1, 'K'*ins_len)
            hgvs_list.append(hgvs_string)
        elif seq_ont == 'ID':
            del_len = len(row.ix['Reference base']) / 3
            hgvs_string = 'p.K{0}_K{1}del{2}'.format(pos, pos+del_len-1, 'K'*del_len)
            hgvs_list.append(hgvs_string)
        elif seq_ont == 'FI' or seq_ont == 'FD':
            hgvs_string = 'p.K{0}fs*'.format(pos)
            hgvs_list.append(hgvs_string)
        elif seq_ont is np.nan:
            hgvs_list.append('p.?')
    return hgvs_list


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--maf',
                        type=str, action='store',
                        help='Filtered MAF file from broad website')
    parser.add_argument('-n', '--non-coding',
                        type=str, default=None,
                        help='Non-coding portion of MAF file from broad'
                        ' website. Only needed if want non-coding.')
    parser.add_argument('-c', '--cravat',
                        type=str, action='store',
                        help='Cravat output which includes amino '
                        'acid info (Variant_Analysis.tsv)')
    parser.add_argument('-o', '--output',
                        type=str, action='store',
                        help='Modified MAF format acceptable for input')
    parser.add_argument('-no', '--non-coding-output',
                        type=str, default=None,
                        help='Modified MAF format for non-coding output')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read in data
    cravat_df = pd.read_csv(opts['cravat'], sep='\t')
    broad_df = pd.read_csv(opts['maf'], sep='\t')
    prev_len = len(broad_df)

    # read non-coding
    if opts['non_coding'] is not None:
        non_coding_df = pd.read_csv(opts['non_coding'], sep='\t')

    # get hgvs strings from cravat output
    hgvs_list = generate_hgvs_syntax(cravat_df)
    broad_df['Protein_Change'] = hgvs_list

    # filter out variants with mappability warning
    broad_df = broad_df[cravat_df['Mappability Warning'].isnull()]
    after_len = len(broad_df)
    print('Before mappability filtering: {0} lines'.format(prev_len))
    print('After mappability filtering: {0} lines'.format(after_len))
    print('Line difference: {0}'.format(prev_len-after_len))

    # rename headers
    rename_cols = {'ttype': 'Tumor_Type',
                   'patient': 'Tumor_Sample',
                   'gene': 'Gene_Symbol',
                   'type': 'Variant_Classification',
                   'chr': 'Chromosome',
                   'ref_allele': 'Reference_Allele',
                   'newbase': 'Tumor_Allele',
                   'pos': 'Start_Position'}
    broad_df = broad_df.rename(columns=rename_cols)

    # add end position column
    broad_df['End_Position'] = broad_df['Start_Position'] + broad_df['Tumor_Allele'].apply(lambda x: len(x) - 1)

    # fix variant classigfication column
    broad_df['Variant_Classification'] = broad_df['Variant_Classification'].apply(fix_variant_type)

    # add chr to chromosome names
    broad_df['Chromosome'] = broad_df['Chromosome'].astype(str).apply(lambda x: 'chr' + x)

    # fix tumor names
    broad_df['Tumor_Type'] = broad_df['Tumor_Type'].apply(fix_tumor_type)
    broad_df['Tumor_Sample'] = broad_df['Tumor_Sample'].apply(fix_tumor_sample)

    cols_of_interest = ['Gene_Symbol', 'Tumor_Sample', 'Tumor_Type',
                        'Chromosome', 'Start_Position',
                        'End_Position', 'Variant_Classification',
                        'Reference_Allele', 'Tumor_Allele',
                        'Protein_Change']

    if opts['non_coding'] is not None:
        non_coding_df = non_coding_df.rename(columns=rename_cols)

        # add end position column
        non_coding_df['End_Position'] = non_coding_df['Start_Position'] + non_coding_df['Tumor_Allele'].apply(lambda x: len(x) - 1)

        # fix variant classigfication column
        non_coding_df['Variant_Classification'] = non_coding_df['Variant_Classification'].apply(fix_variant_type)

        # add chr to chromosome names
        non_coding_df['Chromosome'] = non_coding_df['Chromosome'].astype(str).apply(lambda x: 'chr' + x)

        # fix tumor names
        non_coding_df['Tumor_Type'] = non_coding_df['Tumor_Type'].apply(fix_tumor_type)
        non_coding_df['Tumor_Sample'] = non_coding_df['Tumor_Sample'].apply(fix_tumor_sample)

        # save non-coding
        non_coding_df[cols_of_interest[:-1]].to_csv(opts['non_coding_output'], sep='\t', index=False)

    # output results to file
    broad_df[cols_of_interest].to_csv(opts['output'], sep='\t', index=False)


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
