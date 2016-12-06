Quick Start
===========

This quick start is only meant as a test to check whether 20/20+ has been **properly installed**.
Please see the :ref:`tut-ref` for a more detailed example on the full pipeline that you can modify to use for your own data.

Creating Features
-----------------

First, download the intermediate data files since the output from 
the `probabilistic2020 <http://probabilistic2020.readthedocs.io/en/latest/>`_ package has already been prepared for you.

.. code-block:: bash

    $ wget http://karchinlab.org/data/2020+/pancan_example.tar.gz
    $ tar xvzf pancan_example.tar.gz
    $ cd pancan_example

You do not have to concern yourself with
how these files were generated for the purpose of this quick start.
You should see, however, oncogene.txt for oncogene related features, tsg.txt for tumor suppressor related features, and the summary file named summary_pancan.txt. 

To create the features, use the **features** sub-command for the
2020plus.py script.

.. code-block:: bash

   $ python `which 2020plus.py` features \
        -og-test oncogene.txt \
        -tsg-test tsg.txt \
        --summary summary_pancan.txt \
        -o features_pancan.txt

Prediction
----------

The command for prediction in this quick start example is for pan-cancer data encompassing many cancer types and samples.  
The **classify** sub-command performs the 20/20+ predictions of driver genes.
This step needs the 20/20+ features file already created (features_pancan.txt), and the emprirical 
null distribution file to additionally report p-values/FDR. If an
empirical null distribution file is not provided then only the random
forest scores for prediction. The **empricial null distribution** 
relates the classifier score to a p-value. An example null distribution
**specific** to this pan-cancer dataset is simulated_null_dist.txt.


.. code-block:: bash

   $ python `which 2020plus.py` --out-dir=result_compare classify \
        -f features_pancan.txt \
        -nd simulated_null_dist.txt 

You should see the results in the result_compare/result/r_random_forest_prediction.txt file. It should be the same as the output already provided in result/ directory. Particularly, you should shoulde see 106 TSG scores, 64 oncogene scores, and 197 driver scores (208 unique genes) as significant at a Benjamini-Hochberg FDR of .1.
