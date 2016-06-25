# 20/20+

## About

Next-generation DNA sequencing of the exome has detected hundreds of thousands of small somatic variants (SSV) in cancer. However, distinguishing genes containing driving mutations rather than simply passenger SSVs from a cohort sequenced cancer samples requires sophisticated computational approaches.
20/20+ integrates many features indicative of positive selection to predict oncogenes and tumor suppressor genes from small somatic variants. 
The features capture mutational clustering, conservation, mutation *in silico* pathogenicity scores, mutation consequence types, protein interaction network connectivity, and other covariates (e.g. replication timing).
Contrary to methods based on mutation rate, 20/20+ uses ratio-metric features of mutations by normalizing for the total number of mutations in a gene. This decouples the genes from gene-level differences in background mutation rate. 20/20+ uses monte carlo simulations to evaluate the significance of random forest scores based on an estimated p-value from an empirical null distribution.

## Documentation

[![Documentation Status](http://readthedocs.org/projects/2020plus/badge/?version=latest)](http://2020plus.readthedocs.io/en/latest/?badge=latest)

Please see the [documentation](http://2020plus.readthedocs.io/) on readthedocs.

## Releases

You can download [releases](https://github.com/KarchinLab/2020plus/releases) on github.

* 2020plus v1.0.0 - 5/1/2016 - Initial release

## Installation

[![Build Status](https://travis-ci.org/KarchinLab/2020plus.svg?branch=master)](https://travis-ci.org/KarchinLab/2020plus)

Because 20/20+ internally uses the random forest package in R, you will both need [R](https://www.r-project.org/) and the randomForest library installed. Once R is installed, you can install the random forest package:

```R
> install.packages("randomForest")
```

If you do not have permission to install `randomForest` on the system wide R, you can install in your local user directory by creating an `~/.Renviron` file as the following:

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
* probabilistic2020

To install these packages via `pip` you can use the following command:

```bash
$ pip install -r requirements.txt
```

If you want the exact version 20/20+ was tested on use the `requirements_dev.txt` file and python 2.7. The [probabilistic2020](https://github.com/KarchinLab/probabilistic2020) python package is used to generate features for 20/20+ from mutations in MAF format.
