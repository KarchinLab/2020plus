.. _download-ref:

Download
========

20/20+ releases
---------------

* `2020plus v1.1.0 <https://github.com/KarchinLab/2020plus/archive/v1.1.0.tar.gz>`_ - 11/14/2016 - Improved training procedure and added p-value diagnostic plots
* `2020plus v1.0.3 <https://github.com/KarchinLab/2020plus/archive/v1.0.3.tar.gz>`_ - 10/12/2016 - Fixed error in logging
* `2020plus v1.0.2 <https://github.com/KarchinLab/2020plus/archive/v1.0.2.tar.gz>`_ - 10/03/2016 - Fixed python3 conversion bug
* `2020plus v1.0.1 <https://github.com/KarchinLab/2020plus/archive/v1.0.1.tar.gz>`_ - 6/26/2016 - Added ability to run 20/20+ as a pipeline
* `2020plus v1.0.0 <https://github.com/KarchinLab/2020plus/archive/v1.0.0.tar.gz>`_ - 5/1/2016 - Initial release

Pre-trained classifier
----------------------

We have trained a 20/20+ classifier on pan-cancer data. This can be used to predict on cancer type specific mutations.

Current trained classifier (>= v1.1.0):
* `2020plus_10k.Rdata <http://karchinlab.org/data/2020+/2020plus_10k.Rdata>`_ (NUMSIMULATIONS=10,000, default)
* `2020plus_100k.Rdata <http://karchinlab.org/data/2020+/2020plus_100k.Rdata>`_ (NUMSIMULATIONS=100,000)

Trained classifier for versions 1.0.0-1.0.3 (old):
* `2020plus.Rdata <http://karchinlab.org/data/2020+/2020plus.Rdata>`_

Pan-cancer mutation data
------------------------

* full `pan-cancer <http://karchinlab.org/data/Protocol/pancan-mutation-set-from-Tokheim-2016.txt.gz>`_ data set from:

Tokheim C, Papdopoulis N, Kinzler KW, Vogelstein B, Karchin R (2016) Evaluating the evaluation of cancer driver genes. Submitted `[bioRxiv preprint] <http://biorxiv.org/content/early/2016/06/23/060426>`_

Details about how the mutations were filtered is available `here <https://github.com/KarchinLab/2020plus/blob/master/data/README.rst>`_.

Example data
------------

* Example `pan-cancer <http://karchinlab.org/data/2020+/pancan_example.tar.gz>`_ data set
