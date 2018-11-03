FAQ
===

**Who should I contact if I encounter a problem?**

If you believe your problem may be encountered by other users,
please post the question on `biostars <https://www.biostars.org/>`_.
Check to make sure your question has not been already answered 
by looking at posts with the tag `2020+ <https://www.biostars.org/t/2020+>`_.
Otherwise, create a new post with the 2020+ tag. We will be checking
biostars for questions. You may also contact me directly at
ctokhei1 AT alumni dot jh dot edu.

**Can I use my own custom list of oncogenes/tumor suppressor genes for training?**

Yes, you can change which genes are used for the training list of well-supported oncogenes and tumor suppressor genes. All you need is to change the gene names found
in **data/gene_lists/oncogenes.txt** and **data/gene_lists/tsgs.txt**. This should
only be done for advanced users as the training list of oncogenes and tumor suppressor
genes were established from cancer experts.

**How can I speed up the run time of 20/20+?** 

You can substantially speed up run time by reducing the number of simulations.
This can be done by reducing the NUMSIMULATIONS variable (e.g. from 100000 to 10000) in the `config.yaml` file or specification in the command line of snakemake via `--config NUMSIMULATIONS=10000`. This might result in a slight decrease in prediction performance but may be waranted for large data.

**What happens if silent mutations were not recorded in my data set?** 

Ocassionally, in the literature, studies may only report non-silent mutations
from a sequencing study. If not accounted for, this may bias estimates of 
statistical significance. To make an adjustment for this problem, provide the `drop_silent` option
via the command line: `--config drop_silent="yes"`.
