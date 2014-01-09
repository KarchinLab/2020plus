.. Oncogene Classifier documentation master file, created by
   sphinx-quickstart on Sun Dec  1 16:01:06 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Oncogene Classifier's documentation!
===============================================

This project uses the COSMIC database of mutations to predict which
genes are oncogenes or tumor suppressors. You may want to read the following
papers to understand the code:

* `Cancer Genome Landscape <http://www.sciencemag.org/content/339/6127/1546>`_
* `MutSigCV paper <http://www.nature.com/nature/journal/v499/n7457/full/nature12213.html>`_

It should not be difficult to run my code. The entire pipeline from COSMIC raw
data to predicted genes should only take ~20 minutes.

When I refer to the "20/20 rule" or "vogelstein...", I am referring to the 
"Cancer Genome Landscape" paper methodology. Data about background mutation rates
and olfactory genes come from the MutSigCV paper.


Contents:

.. toctree::
   :maxdepth: 3

   user_doc
   developer_doc
   faq

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

