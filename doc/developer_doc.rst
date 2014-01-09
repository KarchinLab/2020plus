Developer Documentation
=======================

Overview
--------

There are already over 4000 lines of code, as of writting.
However, I have distilled most of the usage to just a handful
of sub-commands of the run.py script in the project root directory.
There are four run.py sub-commands:

1. savedb
2. features
3. data_analysis
4. classify

Python code for each sub-command is located in a self-titled directory.
Lets take the *data_analysis* sub-command as an example:

* data_analysis/

  * python/

    * stats.py
    * . . .

  * results/
  * plots/

Python code is located in the python sub-directory. If figures are generated
then they are stored in the plots/ directory. While text files are stored in 
the results/ directory.

In general, most (if not all) of the file paths for text files are 
found in config files in the *config* directory.

Package dependencies
--------------------

You can find the exact version of every package I used from the 
"requirements.txt" file in the project root directory. Although
you are likely familiar with packages like NumPy or Scipy, you
likely may not know *Pandas* or *Scikit-learn*. A hidden dependency
is the randomForest library in R (I use the rpy2 wrapper).

`Pandas <http://pandas.pydata.org/>`_
+++++++++++++++++++++++++++++++++++++

Pandas is a package that provides an R-like data frame for python.
Learning a little about pandas from their website will help tremendously
in understanding my code. First of all, all I/O operations are handled
through pandas including querying/writing to databases. The following is 
an example of querying sqlite3 and getting a data frame object.

.. code-block:: python

   oncogenes = _utils.oncogene_list  # list of oncogenes
   sql = "SELECT * FROM nucleotide WHERE Gene in " + str(oncogenes)
   df = psql.frame_query(sql, con=conn)  # execute query
   # now do some work . . .

The convenient thing about pandas is that if you can do something to a
data frame in R you likely can accomplish the same thing in pandas.

`scikit-learn <http://scikit-learn.org/stable/>`_
+++++++++++++++++++++++++++++++++++++++++++++++++

Scikit-learn provides the machine learning algorithms in python. Additionally
it has methods for calculating classifier performance and doing cross-validation.

Abbreviations
-------------

I use a considerable amount of abbreviations for names in my code.
The following is a brief, non-inclusive, list that may help when
reading my code:

* **tsg** = tumor suppressor gene
* **onco** = oncogene
* **ct** = count
* **2020** = 20/20 rule

:doc:`code`
-----------

Descriptions of the code can be found on the code description page.
