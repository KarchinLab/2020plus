import pandas as pd
import pandas.rpy.common as com
import rpy2.robjects as ro
import argparse


def venn_diagram(first, second,
                 name1, name2, save_path):
    """Wrapper arround R's VennDiagram."""
    # define function in R to make venn diagram
    ro.r('''venn_diag <- function(df1, df2, save_path){
            library(VennDiagram)
            png()
            venn.diagram(
                x = list(
                    %s = df1$X0,
                    %s = df2$X0
                    ),
                filename = save_path,
                lwd = 4,
                fill = c("cornflowerblue", "darkorchid1"),
                alpha = 0.75,
                label.col = "white",
                cex = 4,
                fontfamily = "serif",
                fontface = "bold",
                cat.col = c("cornflowerblue", "darkorchid1"),
                cat.cex = 3,
                cat.fontfamily = "serif",
                cat.fontface = "bold",
                cat.dist = c(0.03, 0.03),
                cat.pos = c(-20, 14)
                );
            dev.off()
    }''' % (name1, name2))
    venn_diag = ro.r['venn_diag']  # venn diagram function

    # convert to R data frame
    first_rdf = com.convert_to_r_dataframe(first)
    second_rdf = com.convert_to_r_dataframe(second)

    venn_diag(first_rdf, second_rdf, save_path)


def parse_arguments():
    description = 'Perform set operations or a venn diagram'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-f', '--first',
                        type=str, required=True,
                        help='first list')
    parser.add_argument('-s', '--second',
                        type=str, required=True,
                        help='second list')
    parser.add_argument('-nf', '--name-first',
                        type=str, required=True,
                        help='name of first list')
    parser.add_argument('-ns', '--name-second',
                        type=str, required=True,
                        help='name of second list')
    parser.add_argument('-v', '--venn-diagram',
                        type=str, default='',
                        help='Path to save venn diagram (png)')
    args = parser.parse_args()
    return vars(args)


def main(opts):
    # read data
    first = pd.read_csv(opts['fist'], header=None)
    second = pd.read_csv(opts['second'], header=None)

    # make them sets
    first_set = set(first[0].tolist())
    second_set = set(second[0].tolist())

    # print the exact results
    print 'intersect:', len(first_set & second_set)
    print 'first only:', len(first_set - second_set)
    print 'second only:', len(second_set - first_set)
    print 'INTERSECTION'
    print '-' * 20
    for element in (first_set & second_set):
        print element
    print opts['first_name'] + ' only'
    print '-' * 20
    for element in (first_set - second_set):
        print element
    print opts['second_name'] + ' only'
    print '-' * 20
    for element in (second_set - first_set):
        print element

    # plot venn diagram
    venn_diagram(first, second,
                 opts['name_first'], opts['name_second'],
                 opts['venn_diagram'])


if __name__ == "__main__":
    opts = parse_arguments()
    main(opts)
