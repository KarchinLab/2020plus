.. _adv-tut-ref:

Advanced Tutorial
=================

If you use snakemake, you do not need to cover this section.
This advanced tutorial covers all of the individual commands
needed to run the 20/20+ pipeline in a fine grained fashion.
This will make the process more complicated then the regular 
:ref:`tut-ref`.

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
