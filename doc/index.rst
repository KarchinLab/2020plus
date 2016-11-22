.. 20/20+ documentation master file, created by
   sphinx-quickstart on Sun Dec  1 16:01:06 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

20/20+: Ratiometric prediction of cancer driver genes
======================================================

:Author: Collin Tokheim
:Contact: ctokheim AT jhu.edu
:Source code: `GitHub <https://github.com/KarchinLab/2020plus>`_
:Q&A: `Biostars (tag: 2020+) <https://www.biostars.org/t/2020+/>`_ 

Next-generation DNA sequencing of the exome has detected hundreds of thousands of small somatic variants (SSV) in cancer. However, distinguishing genes containing driving mutations rather than simply passenger SSVs from a cohort sequenced cancer samples requires sophisticated computational approaches.
20/20+ integrates many features indicative of positive selection to predict oncogenes and tumor suppressor genes from small somatic variants. 
The features capture mutational clustering, conservation, mutation *in silico* pathogenicity scores, mutation consequence types, protein interaction network connectivity, and other covariates (e.g. replication timing).
Contrary to methods based on mutation rate, 20/20+ uses ratiometric features of mutations by normalizing for the total number of mutations in a gene. This decouples the genes from gene-level differences in background mutation rate. 

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

* `2020plus v1.1.0 <https://github.com/KarchinLab/2020plus/archive/v1.1.0.tar.gz>`_ - 11/21/2016 - Improved training procedure and added p-value diagnostic plots
* `2020plus v1.0.3 <https://github.com/KarchinLab/2020plus/archive/v1.0.3.tar.gz>`_ - 10/12/2016 - Fixed error in logging
* `2020plus v1.0.2 <https://github.com/KarchinLab/2020plus/archive/v1.0.2.tar.gz>`_ - 10/03/2016 - Fixed python3 conversion bug
* `2020plus v1.0.1 <https://github.com/KarchinLab/2020plus/archive/v1.0.1.tar.gz>`_ - 6/26/2016 - Added ability to run 20/20+ as a pipeline
* `2020plus v1.0.0 <https://github.com/KarchinLab/2020plus/archive/v1.0.0.tar.gz>`_ - 5/1/2016 - Initial release

Citation
--------

Tokheim C, Papdopoulis N, Kinzler KW, Vogelstein B, Karchin R (2016) Evaluating the evaluation of cancer driver genes. Submitted `[bioRxiv preprint] <http://biorxiv.org/content/early/2016/06/23/060426>`_
