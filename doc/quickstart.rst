Quick Start
===========

This provides a quick start to running 20/20+ with
the minimum number of steps to execute the statistical test.
It serves as an introction to the workflow for running 20/20+,
with a more detailed tutorial `here <>`_.

Creating Features
-----------------

Creating the features for 20/20+ combines output from
the probabilistic2020 framework. First, download
the data files.

.. code-block:: bash

    $ wget /path/to/og_test
    $ wget /path/to/tsg_test
    $ wget /path/to/summary_file

Where XXX is the probabilistic2020 oncogene output, XXX is the
output of probabilistic2020 tsg, and the summary file
describes many non p-value features.

To create the features, the **features** sub-command for the
2020plus.py script is used.

.. code-block:: bash

   $ python 2020plus.py features \
        -og-test og_file \
        -tsg-test tsg_file \
        --summary summary_file \
        -o features_pancan.txt

Prediction
----------

The first step is to obtain the empirical null distribution for
the random forest scores. You can obtain an already computed
file for this example.

.. code-block:: bash

   $ wget /path/to/nulldist


.. code-block:: bash

   $ python 2020plus.py --out-dir=output_dir classify \
        -f features_pancan.txt \
        -nd nulldist.txt \
