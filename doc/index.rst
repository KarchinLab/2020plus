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

* `2020plus v1.2.0 <https://github.com/KarchinLab/2020plus/archive/v1.2.0.tar.gz>`_ - 3/21/2018 - Change to null distribution simulation
* `2020plus v1.1.3 <https://github.com/KarchinLab/2020plus/archive/v1.1.3.tar.gz>`_ - 8/17/2017 - Bug fixes for different versions of rpy2
* `2020plus v1.1.2 <https://github.com/KarchinLab/2020plus/archive/v1.1.2.tar.gz>`_ - 7/3/2017 - Further bug fixes for latest versions of 20/20+ dependencies
* `2020plus v1.1.1 <https://github.com/KarchinLab/2020plus/archive/v1.1.1.tar.gz>`_ - 5/22/2017 - Bug fixes to work with newest version of pandas
* `2020plus v1.1.0 <https://github.com/KarchinLab/2020plus/archive/v1.1.0.tar.gz>`_ - 11/21/2016 - Improved training procedure and added p-value diagnostic plots
* `2020plus v1.0.3 <https://github.com/KarchinLab/2020plus/archive/v1.0.3.tar.gz>`_ - 10/12/2016 - Fixed error in logging
* `2020plus v1.0.2 <https://github.com/KarchinLab/2020plus/archive/v1.0.2.tar.gz>`_ - 10/03/2016 - Fixed python3 conversion bug
* `2020plus v1.0.1 <https://github.com/KarchinLab/2020plus/archive/v1.0.1.tar.gz>`_ - 6/26/2016 - Added ability to run 20/20+ as a pipeline
* `2020plus v1.0.0 <https://github.com/KarchinLab/2020plus/archive/v1.0.0.tar.gz>`_ - 5/1/2016 - Initial release

Citation
--------

Collin J. Tokheim, Nickolas Papadopoulos, Kenneth W. Kinzler, Bert Vogelstein, and Rachel Karchin. Evaluating the evaluation of cancer driver genes. PNAS 2016 ; published ahead of print November 22, 2016, doi:10.1073/pnas.1616440113
