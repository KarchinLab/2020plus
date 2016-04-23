20/20+ Classifier
-----------------

20/20+ integrates many features indicative of positive selection to predict oncogenes and tumor suppressor genes from small somatic variants. 
The features capture mutational clustering, conservation, mutation *in silico* pathogenicity scores, mutation consequence types, protein interaction network connectivity, and other covariates (e.g. replication timing).
Contrary to methods based on mutation rate, 20/20+ uses ratio-metric features of mutations by normalizing for the total number of mutations in a gene. This decouples the genes from gene-level differences in background mutation rate. 20/20+ uses monte carlo simulations to evaluate the significance of random forest scores based on an estimated p-value from an empirical null distribution.


Installation
------------

[![Build Status](https://travis-ci.com/ctokheim/2020plusClassifier.svg?token=KhnctpTdxNuuZ9Z1kcsg&branch=master)](https://travis-ci.com/ctokheim/2020plusClassifier)

Because 20/20+ internally uses the random forest package in R, you will both need [R](https://www.r-project.org/) and the randomForest library installed. Once R is installed, you can install the random forest package:

```R
> install.packages("randomForest")
```

If you do not have permission to install `randomFroest` on the system wide R, you can install in your local user directory by creating an `~/.Renviron` file as the following:

```
R_LIBS=~/Rlibs
```

Where, in this case, the R libraries will be installed in the `~/Rlibs` directory.

20/20+ also requires the following python packages:

* numpy
* scipy
* pandas>=0.17.0
* scikit-learn
* rpy2

To install these packages via `pip` you can use the following command:

```bash
$ pip install -r requirements.txt
```

If you want the exact version 20/20+ was tested on use the `requirements_dev.txt` file and python 2.7. For generating your own input features for 20/20+ please install the [probabilistic2020](https://github.com/KarchinLab/probabilistic2020) python package.
