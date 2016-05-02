Quick Start
===========

This provides a quick start to running 20/20+ with
the minimum number of steps to execute the statistical test.
It serves as an introction to the workflow for running 20/20+,
see the :ref:`tut-ref` for a more detailed example.

Creating Features
-----------------

Creating the features for 20/20+ combines output from
the  `probabilistic2020 <http://probabilistic2020.readthedocs.org/>`_
package. First, download
the data files.

.. code-block:: bash

    $ wget http://karchinlab.org/data/2020+/pancan_example.tar.gz
    $ tar xvzf pancan_example.tar.gz
    $ cd pancan_example

Where oncogene.txt is probabilistic2020 oncogene output, tsg.txt is the
output of probabilistic2020 tsg, and the summary file named summary_pancan.txt
describes many non p-value features.

To create the features, use the **features** sub-command for the
2020plus.py script.

.. code-block:: bash

   $ python 2020plus.py features \
        -og-test oncogene.txt \
        -tsg-test tsg.txt \
        --summary summary_pancan.txt \
        -o features_pancan.txt

Prediction
----------

The first step is to obtain the empirical null distribution for
the random forest scores. The **empricial null distribution** 
relates the classifier score to a p-value. An example null distribution
specific to this pan-cancer dataset is simulated_null_dist.txt.

The **classify** sub-command performs the 20/20+ predictions of driver genes.
This step needs the 20/20+ features file already created, and the emprirical 
null distribution file to additionally report p-values/FDR. If an
empirical null distribution file is not provided then only the random
forest scores for prediction.

.. code-block:: bash

   $ python 2020plus.py --out-dir=result_compare classify \
        -f features_pancan.txt \
        -nd simulated_null_dist.txt 

You should see the results in the result_compare/classify/result/r_random_forest_prediction.txt file. It should be the same as the output already provided in result/ directory. Particularly, you should shoulde see 106 TSG scores, 64 oncogene scores, and 197 driver scores (208 unique genes) as significant at a Benjamini-Hochberf FDR of .1.
