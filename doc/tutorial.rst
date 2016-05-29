.. _tut-ref:

Tutorial
========

Background
----------

20/20+ uses the random forest algorithm to predict cancer driver genes.
Because 20/20+ uses supervised machine learning, cancer driver gene prediction depends
on a list of defined oncogenes and tumor suppressor genes used for training. The list is already
pre-defined from the `Cancer Genome Landscapes <http://www.ncbi.nlm.nih.gov/pubmed/23539594>`_ paper. 
10-fold cross-validation is performed internally within 20/20+ to avoid overfitting.
Included in the results are a score specifically for oncogene or tumor suppressor gene,
and a cumulative driver gene score. The score represents the fraction of decision
trees in the random forest which voted for the particular class (oncogene, TSG, driver (either oncogene or TSG), passenger). Scores nearer one indicate stronger evidence for the gene as a cancer driver.

There are two ways to use 20/20+ scores. One is to just obtain cancer driver scores for genes, which allows ranking of genes in a prioritized manner. The second way evaluates the statistical significance of 
the cancer driver score. A p-value and associated Benjamini-Hochberg false discovery rate
will be reported, but this requires establishing the null distribution for random forest scores.
20/20+ uses an empirical null distribution by scoring simulated mutations in genes.
This extra step requires additional work and computational resources.

Creating feature matrix
-----------------------

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
+++++++++++++++++++++++++++++++++

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
++++++++++++++++

A single feature file ("features.txt") containing all three of the above commands is created
by the **2020plus.py features** command.

.. code-block:: bash

   $ python 2020plus.py features \
        -og-test oncogene.txt \
        -tsg-test tsg.txt \
        --summary summary.txt \
        -o features.txt


Predicting cancer driver genes
------------------------------

Scores only
+++++++++++


Statistical significance
++++++++++++++++++++++++

.. image:: /images/final_result.png
    :scale: 50%
    :align: center

Creating null distribution
##########################

The first step is to obtain a trained classifier on the observed data.
You can skip this step if you download an already trained classifier
used in Tokheim et al. (`here <>`_).

.. image:: /images/pancan_trained_classifier.png
    :scale: 50%
    :align: center

Next, create simulated data.

.. image:: /images/simulated_features.png
    :align: center

Finally, score the simulations to obtain an empirical null distribution.

.. image:: /images/null_distribution.png
    :scale: 50%
    :align: center
