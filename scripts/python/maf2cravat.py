import csv
import argparse


def parse_arguments():
    descript = ('Convert MAF files '
                ' into format compatible with cravat')
    parser = argparse.ArgumentParser(description=descript)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-s', '--saturation-paper',
                       action='store_true',
                       help="Parse cancer gene saturation paper")
    group.add_argument('-t', '--tuson-paper',
                       action='store_true',
                       help="Parse tuson paper")
    args, unknown = parser.parse_known_args()
    return vars(args), unknown


def main(opts, argv):
    maf_file, cravat_file = argv[0], argv[1]
    out_list = []
    with open(maf_file) as handle, open(cravat_file, 'wb') as mywriter:
        maf_list = list(csv.reader(handle, delimiter='\t'))
        header = maf_list.pop(0)
        col2ix = {c: i for i, c in enumerate(header)}
        for k, line in enumerate(maf_list):
            if opts['saturation_paper']:
                # add the extra base for cravat if an indel
                is_ins = line[col2ix['classification']] == 'INS' or '-' in line[col2ix['ref_allele']]
                is_del = line[col2ix['classification']] == 'DEL' or '-' in line[col2ix['newbase']]
                is_indel = is_ins or is_del
                if is_indel:
                    ref_base = 'A' + line[col2ix['ref_allele']].replace('-', '')
                    new_base = 'A' + line[col2ix['newbase']].replace('-', '')
                else:
                    ref_base = line[col2ix['ref_allele']]
                    new_base = line[col2ix['newbase']]

                # fix chr names
                mychr = line[col2ix['chr']]
                mychr = 'X' if mychr == '23' else mychr
                mychr = 'Y' if mychr == '24' else mychr
                mychr = 'chr' + mychr

                # get position
                mypos = line[col2ix['pos']]
            else:
                pass

            # append results
            tmp_list = [k, mychr, mypos, '+',
                        ref_base, new_base]
            out_list.append(tmp_list)
        csv.writer(mywriter, delimiter='\t').writerows(out_list)


if __name__=="__main__":
    opts, argv = parse_arguments()
    main(opts, argv)
