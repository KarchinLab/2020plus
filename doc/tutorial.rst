.. _tut-ref:

Tutorial
========

Since 20/20+ uses supervised machine learning, cancer driver gene prediction depends
on a list of defined oncogenes and tumor suppressor genes. The list is already
pre-defined from the Cancer Genome Landscapes paper. Thus what is needed 
is to populate a feature matrix for prediction and create an empirical null 
distribution to define significance.

Creating feature matrix
-----------------------

.. image:: /images/pancan_feat_matrix.png
    :align: center

Creating null distribution
--------------------------

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

Predicting cancer driver genes
------------------------------

.. image:: /images/final_result.png
    :scale: 50%
    :align: center
