Installation
------------

.. image:: https://travis-ci.org/KarchinLab/2020plus.svg?branch=master
    :target: https://travis-ci.org/KarchinLab/2020plus

Releases
~~~~~~~~

20/20+ can be downloaded on `github <https://github.com/KarchinLab/2020plus/releases>`_.

Package requirements
~~~~~~~~~~~~~~~~~~~~

Once you have downloaded the source code for 20/20+, please add the directory to your :code:`PATH` variable.

We recommend that you install the dependencies for 20/20+ through `conda <https://conda.io/miniconda.html>`_. Once conda is installed, setting up the environment is done as follows:

.. code-block:: bash

    $ conda env create -f environment_python.yml  # install dependencies for python
    $ source activate 2020plus  # activate the 20/20+ conda environment
    $ conda install r r-randomForest rpy2  # install the R related dependencies

Every time you wish to run 20/20+, you will then need to activate the "2020plus" conda environment.

.. code-block:: bash

    $ source activate 2020plus

The 20/20+ conda environment can also be deactivated.

.. code-block:: bash

    $ source deactivate 2020plus

Check your PATH variable
~~~~~~~~~~~~~~~~~~~~~~~~

Make sure that you have add the 20/20+ directory to your `PATH` variable. We recommend you add this line to your bashrc file.

.. code-block:: bash

    export PATH=$PATH:/path/to/2020plus

Where "/path/to/2020plus" represents the path where **you** placed 20/20+. If you have done this correctly, the following command should print the location of the 2020plus.py script.

.. code-block:: bash

    $ which 2020plus.py
