.. _download-ref:

Download
========

20/20+ releases
---------------

* `2020plus v1.2.3 <https://github.com/KarchinLab/2020plus/archive/v1.2.3.tar.gz>`_ - 4/6/2019 - Minor change for installation procedure
* `2020plus v1.2.2 <https://github.com/KarchinLab/2020plus/archive/v1.2.2.tar.gz>`_ - 9/10/2018 - Added option to handle mutational data sets where silent mutations are not reported
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

* `Pre-computed scores <https://www.dropbox.com/scl/fi/o6eaih9d3rr9ztms2pa33/scores.tar.gz?rlkey=ekih9qstzfncn8811935a7ghb&st=58ife2e0&dl=1>`_ data set
* `Reference SNVBox transcripts <https://www.dropbox.com/scl/fi/jnucugcu4qslb8vw0s9ry/snvboxGenes.bed?rlkey=yc3gqu4msx0wgqb6wo149rpok&st=dklww6yv&dl=1>`_ in BED format

Pre-trained classifier
----------------------

We have trained a 20/20+ classifier on pan-cancer data. This can be used to predict on cancer type specific mutations.

Current trained classifier (>= v1.1.0):

* `2020plus_10k.Rdata <https://www.dropbox.com/scl/fi/zv3twjeii2ghxtgy4f9o5/2020plus_10k.Rdata?rlkey=yu8i09tuuf6bbgfzsm7dcynp8&st=yhr9091b&dl=1>`_ (NUMSIMULATIONS=10,000, default)
* `2020plus_100k.Rdata <https://www.dropbox.com/scl/fi/kouu7021tn2zr7l55ws10/2020plus_100k.Rdata?rlkey=ravsamve3zwfa642yxdifm5i3&st=tigmaqba&dl=1>`_ (NUMSIMULATIONS=100,000)

Trained classifier for versions 1.0.0-1.0.3 (old):

* `2020plus.Rdata <https://www.dropbox.com/scl/fi/oaspjxao9fm70wx9tmscm/2020plus.Rdata?rlkey=n82vxs1rak3y81quenoqxbaz7&st=lfdt2vwh&dl=1>`_

Pan-cancer mutation data
------------------------

* full `pan-cancer <https://www.dropbox.com/scl/fi/8ob367fu9ztplyx4mmcj0/pancan-mutation-set-from-Tokheim-2016.txt.gz?rlkey=nxxwkotnuggw2ptinjbbfvp96&st=lw6ah2ip&dl=1>`_ data set from:

Collin J. Tokheim, Nickolas Papadopoulos, Kenneth W. Kinzler, Bert Vogelstein, and Rachel Karchin. Evaluating the evaluation of cancer driver genes. PNAS 2016 ; published ahead of print November 22, 2016, doi:10.1073/pnas.1616440113

Details about how the mutations were filtered is available `here <https://github.com/KarchinLab/2020plus/blob/master/data/README.rst>`_.

Example data
------------

* Example `pan-cancer <https://www.dropbox.com/scl/fi/phv6sidw55qwiqq0g551a/pancan_example.tar.gz?rlkey=psnmy08lcssjj60gagowhwrf8&st=u2ndd173&dl=1>`_ data set
