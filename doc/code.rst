Code
====

The code is written in almost 100% python. This has a clear advantage of
preventing usage of "hacked" scripts. In fact, the code is meant to
be used through a single python script named run.py. If a module has
a main method then generally that is the only function that is called
externally from a different module. The one exception is the use of 
RPy2 to run R's random forest algorith implemented in the randomForest
library.

:doc:`data_analysis`
--------------------

The data_analysis folder contain python modules that calculate statistics
of the COSMIC database. The data_analysis also contains code for plotting
the resulting data. Code for generating the feature matrix for classification
is also located within the data_analysis folder.

:doc:`classify`
---------------

The classify folder contains python modules implementing classification method.
This includes:

* 20/20 rule
* random forest
* naive bayes

:doc:`utils`
------------

The utils folder contains a wide variety of python modules that are generally useful.
This includes the basic plot.py module which does the plotting in matplotlib. In addition,modules like amino_acid and nucleotide parse the HGVS mutation syntax found in the 
COSMIC database.

tests
-----

The tests folder contains python modules for unit tests. All the modules start with
and follow the module/function naming rules for the nose unit testing package.
To run unit tests just run "nosetests tests/module.py".
