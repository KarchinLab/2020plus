.. _tut-ref:

Tutorial
========

Technical background
--------------------

20/20+ uses the random forest algorithm to predict cancer driver genes.
Because 20/20+ uses supervised machine learning, cancer driver gene prediction depends
on a list of defined oncogenes and tumor suppressor genes used for training. The list is already
pre-defined from the `Cancer Genome Landscapes <http://www.ncbi.nlm.nih.gov/pubmed/23539594>`_ paper. 
Included in the results are a score specifically for oncogene or tumor suppressor gene,
and a cumulative driver gene score. The score represents the fraction of decision
trees in the random forest which voted for the particular class (oncogene, TSG, driver (either oncogene or TSG), passenger). Scores nearer one indicate stronger evidence for the gene as a cancer driver.  There are two ways to use 20/20+ scores. One is to just obtain cancer driver scores for genes, which allows ranking of genes in a prioritized manner. The second way evaluates the statistical significance of 
the cancer driver score. A p-value and associated Benjamini-Hochberg false discovery rate
will be reported, but this requires establishing the null distribution for random forest scores.
20/20+ uses an empirical null distribution by scoring simulated mutations in genes.
This extra step requires additional work and computational resources.
2/20+ can be applied to pan-cancer and tumor type specific data. 10-fold cross-validation is performed internally within 20/20+ to avoid overfitting.

20/20+ pipeline
---------------

The easiest way to run the entire 20/20+ pipeline from somatic mutations to cancer
driver gene prediction is to use `snakemake <https://bitbucket.org/snakemake/snakemake/wiki/Home>`_. We have created a **Snakefile** that will run the multiple steps needed to get the final results. You first will need to install snakemake (requires python 3.X), so please see the snakemake `installation instructions <https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-installation>`_. If instead you are using python 2.7, please see the section on :ref:`adv-tut-ref`.

There are two ways to perform predictions with 20/20+. Either in a pan-cancer setting where mutations from several cancer types are aggregated together, or predicting cancer type specific driver genes by using a 20/20+ model previously trained on pan-cancer data.

Pan-cancer analysis
+++++++++++++++++++

.. note:: The pan-cancer tutorial is more computationally intensive and the 
          run time will take a while even on a computer cluster. However, it does 
          demonstrate the correct usage of the 20/20+ pipeline, which is quicker for
          cancer type specific data sets.

The **snakemake -s Snakefile predict** command will perform predictions on pan-cancer
data. Here, it is assumed you are in the 2020plus directory where the Snakefile is located.

.. _prep-data-ref:

Preparing data
##############

To set up this tutorial example you will first need to download the data.
This step is needed regardless of whether your doing predictions on pan-cancer
or cancer type specific data.

.. code-block:: bash

   $ mkdir data
   $ cd data
   $ wget http://karchinlab.org/data/Protocol/pancan-mutation-set-from-Tokheim-2016.txt.gz  # download mutations
   $ wget http://karchinlab.org/data/2020+/snvboxGenes.bed  # download transcript annotation
   $ wget http://karchinlab.org/data/2020+/scores.tar.gz  # download pre-computed scores
   $ gunzip pancan-mutation-set-from-Tokheim-2016.txt.gz 
   $ mv pancan-mutation-set-from-Tokheim-2016.txt mutations.txt  # rename file
   $ tar xvzf scores.tar.gz
   $ cd ..

The somatic mutations in the mutations.txt file is a MAF-like format described `here <http://probabilistic2020.readthedocs.io/en/latest/tutorial.html#mutations>`_. snvboxGenes.bed
contains reference transcripts for each gene (described `here <http://probabilistic2020.readthedocs.io/en/latest/tutorial.html#gene-bed-file>`_), and the scores directory contains pre-computed scores (described `here <http://probabilistic2020.readthedocs.io/en/latest/tutorial.html#pre-computed-scores-optional>`_). You will also need to create a FASTA file
of gene sequences by the following `these instructions <http://probabilistic2020.readthedocs.io/en/latest/tutorial.html#gene-fasta>`_.

Running 20/20+
##############

.. note:: You can substantially speed up run time by reducing the number of simulations.
          This can be done by reducing the NUMSIMULATIONS variable (e.g. from 100000 to 10000) in the `config.yaml` file or specification in the command line of snakemake via `--config NUMSIMULATIONS=10000`. This might result in a slight decrease in prediction performance but may be waranted for large data.

By default, the data is assumed to be located in the "data/" directory and mutations are
"data/mutations.txt". You can change the default by editing the config.yaml file.
However you can also override the default from the command line by specifying
variables with the **--config** argument. The following command executes
the 20/20+ on a local machine.

.. code-block:: bash

   $ snakemake -s Snakefile predict -p --cores 1 \
        --config mutations="data/mutations.txt" output_dir="output"

The **--cores** argument specifies the number of computer cores that are allowable
to be used at a given time.

It is generally recommended to run 20/20+ on a cluster to parallelize
calculations. The below command will execute
the 20/20+ pipeline on an SGE computer cluster using qsub.

.. code-block:: bash

   $ snakemake -s Snakefile predict -p -j 999 -w 10 --max-jobs-per-second 1 \
        --config mutations="data/mutations.txt" output_dir="output" \
        --cluster-config cluster.yaml \
        --cluster "qsub -cwd -pe smp {threads} -l mem_free={cluster.mem},h_vmem={cluster.vmem} -v PATH=$PATH"

In this example, the output will be saved in the "output" directory as specified by the
output_dir parameter (also changeable in config.yaml). The **--cluster** argument
specifies the command prefix for submitting to your cluster job scheduler.
In the above example, **qsub** is used for the SGE scheduler, but this obviously
is cluster specific and therefore you should look up the manual for your cluster.
Of importance, though, is that certain template values can be inserted in to
the job submission. Templated values are denoted by curly braces, and are used
to set the number of threads ("{threads}") and memory ("{cluster.mem}" and "{cluster.vmem}").
Templated values with "cluster." are specified in the cluster config file (cluster.yaml; **--cluster-config** argument). It is also recommended that your PATH environmental variable
is passed into the cluster job submission so that you do not receive a command not found
error. The "-j" argument can restrict the number of concurrent jobs submitted to the cluster,but in our case we use 999 to let the cluster job scheduler to identify which jobs get executed.
The "-w 10 --max-jobs-per-second 1" parameters are issued to avoid overly quick 
job submissions to the cluster.


20/20+ output
#############

Like in the quick start, you will find the result in output/results/r_random_forest_prediction.txt. There will be a p-value/q-value for the oncogene, tumor suppressor gene, and driver
score. The file will also contain all of the features used for prediction.

Cancer type specific analysis
+++++++++++++++++++++++++++++

When performing predictions on cancer type specific mutations, a pre-trained
20/20+ classifier based on pan-cancer data is used to make predictions. The 
first step is to download the `pre-trained 20/20+ <http://karchinlab.org/data/2020+/2020plus.Rdata>`_. Associated data should be collected like for the pan-cancer :ref:`prep-data-ref` section. Instead of using the **predict** command, the **snakemake -s Snakefile pretrained_predict** command should be used. In the below example command, we use the command for a local machine, but as like in the previous example, it can be adopted to run on a cluster.

.. note:: Care should be taken if you intend to predict on samples which were
          used for training the 20/20+ random forest (e.g. predicting on TCGA data).
          This could result in over-fitting, and may require training a 20/20+ random forest yourself (see :ref:`train-ref`) without the overlapping samples.

.. code-block:: bash

   $ snakemake -s Snakefile pretrained_predict -p --cores 1 \
        --config mutations="data/my_cancer_specific_mutations.txt" output_dir="output" trained_classifier="data/2020plus.Rdata"

The difference with the previous pan-cancer command is that the mutations ("data/my_cancer_specific_mutations.txt") are from a single cancer type, and the pre-trained classifier is specified with the **trained_classifier** option. In this case the pre-trained 20/20+ classifier was assumed to be placed into the data directory.

.. _train-ref:

Train a 20/20+ classifier
+++++++++++++++++++++++++

You can also train your own 20/20+ model to predict on new data (e.g. new cancer type specific data) using the **train** command. Training should be performed on a pan-cancer collection of mutations. This either could be those `mutations <http://karchinlab.org/data/Protocol/pancan-mutation-set-from-Tokheim-2016.txt.gz>`_ used in our evaluation or a new collected set. Note, the provided `pre-trained classifier <http://karchinlab.org/data/2020+/2020plus.Rdata>`_ is already trained on the mutations linked in the previous sentence. The file format for mutations is described `here <http://probabilistic2020.readthedocs.io/en/latest/tutorial.html#mutations>`_. Like above, the command can be easily modified to run on a cluster.

.. code-block:: bash

   $ snakemake -s Snakefile train -p --cores 1 \
        --config mutations="data/my_pancancer_mutations.txt" output_dir="output" 

where "data/my_pancancer_mutations.txt" is the file containing small somatic mutations and the trained 20/20+ model will be saved as "output/2020plus.Rdata".
