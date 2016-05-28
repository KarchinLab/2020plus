.. 20/20+ documentation master file, created by
   sphinx-quickstart on Sun Dec  1 16:01:06 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

20/20+: Ratio-metric prediction of cancer driver genes
======================================================

:Author: Collin Tokheim
:Contact: ctokheim AT jhu.edu
:Source code: `GitHub <https://github.com/KarchinLab/2020plus>`_
:Q&A: `Biostars (tag: 2020+) <https://www.biostars.org/t/2020+/>`_ 

Next-generation DNA sequencing of the exome has detected hundreds of thousands of small somatic variants (SSV) in cancer. However, distinguishing genes containing driving mutations rather than simply passenger SSVs from a cohort sequenced cancer samples requires sophisticated computational approaches.
20/20+ integrates many features indicative of positive selection to predict oncogenes and tumor suppressor genes from small somatic variants. 
The features capture mutational clustering, conservation, mutation *in silico* pathogenicity scores, mutation consequence types, protein interaction network connectivity, and other covariates (e.g. replication timing).
Contrary to methods based on mutation rate, 20/20+ uses ratio-metric features of mutations by normalizing for the total number of mutations in a gene. This decouples the genes from gene-level differences in background mutation rate. 20/20+ uses monte carlo simulations to evaluate the significance of random forest scores based on an estimated p-value from an empirical null distribution.


Contents:

.. toctree::
   :maxdepth: 3

   download
   installation
   quickstart
   tutorial
   faq

Releases
--------

* `2020plus v1.0.0 <https://github.com/KarchinLab/2020plus/archive/v1.0.0.tar.gz>`_ - 5/1/2016 - Initial release
