.. _download-ref:

Download
========

20/20+ releases
---------------

* `2020plus v1.2.1 <https://github.com/KarchinLab/2020plus/archive/v1.2.1.tar.gz>`_ - 8/2/2018 - Fixed bug where configuration file would not load
* `2020plus v1.2.0 <https://github.com/KarchinLab/2020plus/archive/v1.2.0.tar.gz>`_ - 3/21/2018 - Change to null distribution simulation
* `2020plus v1.1.3 <https://github.com/KarchinLab/2020plus/archive/v1.1.3.tar.gz>`_ - 8/17/2017 - Bug fixes for different versions of rpy2
* `2020plus v1.1.2 <https://github.com/KarchinLab/2020plus/archive/v1.1.2.tar.gz>`_ - 7/3/2017 - Further bug fixes for latest versions of 20/20+ dependencies
* `2020plus v1.1.1 <https://github.com/KarchinLab/2020plus/archive/v1.1.1.tar.gz>`_ - 5/22/2017 - Bug fixes to work with newest versions of pandas
* `2020plus v1.1.0 <https://github.com/KarchinLab/2020plus/archive/v1.1.0.tar.gz>`_ - 11/21/2016 - Improved training procedure and added p-value diagnostic plots
* `2020plus v1.0.3 <https://github.com/KarchinLab/2020plus/archive/v1.0.3.tar.gz>`_ - 10/12/2016 - Fixed error in logging
* `2020plus v1.0.2 <https://github.com/KarchinLab/2020plus/archive/v1.0.2.tar.gz>`_ - 10/03/2016 - Fixed python3 conversion bug
* `2020plus v1.0.1 <https://github.com/KarchinLab/2020plus/archive/v1.0.1.tar.gz>`_ - 6/26/2016 - Added ability to run 20/20+ as a pipeline
* `2020plus v1.0.0 <https://github.com/KarchinLab/2020plus/archive/v1.0.0.tar.gz>`_ - 5/1/2016 - Initial release

Necessary data files
--------------------

* `Pre-computed scores <http://karchinlab.org/data/2020+/scores.tar.gz>`_ data set
* `Reference SNVBox transcripts <http://karchinlab.org/data/2020+/snvboxGenes.bed>`_ in BED format

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

Collin J. Tokheim, Nickolas Papadopoulos, Kenneth W. Kinzler, Bert Vogelstein, and Rachel Karchin. Evaluating the evaluation of cancer driver genes. PNAS 2016 ; published ahead of print November 22, 2016, doi:10.1073/pnas.1616440113

Details about how the mutations were filtered is available `here <https://github.com/KarchinLab/2020plus/blob/master/data/README.rst>`_.

Example data
------------

* Example `pan-cancer <http://karchinlab.org/data/2020+/pancan_example.tar.gz>`_ data set
