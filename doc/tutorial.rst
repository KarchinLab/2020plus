.. _tut-ref:

Tutorial
========

Background
----------

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
driver gene prediction is to use `snakemake <https://bitbucket.org/snakemake/snakemake/wiki/Home>`_. We have created a **Snakefile** that will run the multiple steps needed to get the final results. You first will need to install snakemake (requires python 3.X), so please see the snakemake `installation instructions <https://bitbucket.org/snakemake/snakemake/wiki/Documentation#markdown-header-installation>`_. If instead you are using python 2.7, please see the section on :ref:`ind-cmd-ref`.

There are two ways to perform predictions with 20/20+. Either in a pan-cancer setting where mutations from several cancer types are aggregated together, or predicting cancer type specific driver genes by using 20/20+ previously trained on pan-cancer data.

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
score.

Cancer type specific analysis
+++++++++++++++++++++++++++++

When performing predictions on cancer type specific mutations, a pre-trained
20/20+ classifier based on pan-cancer data is used to make predictions. The 
first step is to download the `pre-trained 20/20+ <http://karchinlab.org/data/2020+/2020plus.Rdata>`_. Associated data should be collected like for the pan-cancer :ref:`prep-data-ref` section. Instead of using the **predict** command, the **snakemake -s Snakefile pretrained_predict** command should be used. In the below example command, we use the command for a local machine, but as like in the previous example, it can be adopted to run on a cluster.

.. code-block:: bash

   $ snakemake -s Snakefile pretrained_predict -p --cores 1 \
        --config mutations="data/my_cancer_specific_mutations.txt" output_dir="output" trained_classifier="data/2020plus.Rdata"

The difference with the previous pan-cancer command is that the mutations ("data/my_cancer_specific_mutations.txt") are from a single cancer type, and the pre-trained classifier is specified with the **trained_classifier** option. In this case the pre-trained 20/20+ classifier was assumed to be placed into the data directory.

.. _ind-cmd-ref:

Running individual commands
---------------------------

If you use snakemake, you do not need to cover this section.

Creating feature matrix
+++++++++++++++++++++++

20/20+ uses a total of 24 features encompassing mutational clustering,
functional impact bias, evolutionary conservation, composition of mutation consequence types,
protein interaction network, and mutation rate covariates (e.g. replication timing).
The `probabilistic2020 package <http://probabilistic2020.readthedocs.org>`_ is used to 
generate three files that will then be combined for all of the need features.
In the quick start, these files were already provided for you. A flow chart diagram of the 
steps are shown below. It is recommended that you examine the `probabilistic2020 documentation <http://probabilistic2020.readthedocs.org>`_
before continuing.

.. image:: /images/pancan_feat_matrix.png
    :align: center

Running probabilistic2020 package
#################################

The first step is to compute several summary statistics of the composition of mutations in each
gene. This is done with the **mut_annotate** command with the **--summary** flag as follows.

.. code-block:: bash

   $ mut_annotate --summary \
        -i genes.fa \
        -b genes.bed \
        -s score_dir \
        -m mutations.txt \
        -o summary.txt

Where genes.fa is your gene FASTA file for your reference transcripts in genes.bed, mutations.txt is your MAF file containing mutations, score_dir is the directory containing the pre-computed `VEST <http://www.ncbi.nlm.nih.gov/pubmed/23819870>`_ and evolutionary conservation scores, and summary.txt is the output file containing the features. The pre-computed scores
are based on hg19 and the reference transcript annotation from SNVBox. If you are not
using either, then the parameter may be left empty resulting in no features for VEST or evolutionary conservation.

The next two steps involve running a statistical test for features commonly associated
with oncogenes (mutational clustering and mutation functional impact bias) and
TSGs (high proportion of inactivating mutations). The p-values from these
statistical tests are used as features in the 20/20+ predictions. The first 
command performs a test for TSGs.

.. code-block:: bash

   $ probabilistic2020 tsg \
        -i genes.fa \
        -b genes.bed \
        -m mutations.txt \
        -p 1 \
        -n 100000 \
        -o tsg.txt

Because evaluating statistical significance for large datasets can be computationally 
intensive, there are a couple parameters which can speed up calculations. 
In the above example **-p 1** indicates 1 processes should be used, but can be increased to parallelize the calculations (e.g. **-p 10** to split calculations on 10 processes). 
Since the p-value is obtained by simulations, higher number of simulations (-n parameter) means 
increased precision in the reported p-value, but at the cost of increased run time.
The recommended default is **-n 100000** simulations but can be tweaked to obtain a good
balance in run-time and precision of p-value. Generally it is recommended to run
the command on a server with multiple processes (**-p** parameter) rather than lowering
the number of simulations.

The other command intended for oncogenes is the following.

.. code-block:: bash

   $ probabilistic2020 oncogene \
        -i genes.fa \
        -b genes.bed \
        -m mutations.txt \
        -s score_dir \
        -p 1 \
        -n 100000 \
        -o oncogene.txt

Like the mut_annotate command, providing the directory ("score_dir") for pre-computed VEST and evolutionary
conservation scores is optional. However, this will result in not including some important features
into 20/20+ likely decreasing performance.

Merging features
################

A single feature file ("features.txt") containing all three of the above commands is created
by the **2020plus.py features** command.

.. code-block:: bash

   $ python 2020plus.py features \
        -og-test oncogene.txt \
        -tsg-test tsg.txt \
        --summary summary.txt \
        -o features.txt

Predicting cancer driver genes
++++++++++++++++++++++++++++++

Scores only
###########

If interested in only scoring genes, then the next step
is prediction. This is performed with the **2020plus.py classify**
command.

.. code-block:: bash

   $ python 2020plus.py --out-dir=myresult_dir classify -f features.txt 

Where myresult_dir is the directory where results are saved, and features.txt
is the feature file from the **2020plus.py features** command.

Statistical significance
########################

Obtaining a p-value for driver scores requires creating an empirical null distribution 
for use in the prediction step, as diagrammed below.

.. image:: /images/final_result.png
    :scale: 50%
    :align: center

Creating null distribution
~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step is to obtain a trained classifier on the observed data.
You can skip this step if you download an already trained classifier
used in Tokheim et al. (`here <http://karchinlab.org/data/2020+/2020plus.Rdata>`_). The procedure is diagrammed below, and
is critical that a pan-cancer mutation data set is used for training.

.. image:: /images/pancan_trained_classifier.png
    :scale: 50%
    :align: center

Saving a trained classifier is done using the **2020plus.py train** command.

.. code-block:: bash

   $ python 2020plus.py train -f features.txt -r classifier.Rdata 

Where features.txt is the feature file from pan-cancer mutation data set, and
classifier.Rdata is the trained 20/20+ classifier file.

The next step is to create simulated mutations that mimic the random accumulation
of passenger mutations. A diagram of the steps is shown below.

.. image:: /images/simulated_features.png
    :align: center

Finally, score the simulations to obtain an empirical null distribution.

.. image:: /images/null_distribution.png
    :scale: 50%
    :align: center

Prediction
~~~~~~~~~~
