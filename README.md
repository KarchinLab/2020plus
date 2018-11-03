# 20/20+

## About

Next-generation DNA sequencing of the exome has detected hundreds of thousands of small somatic variants (SSV) in cancer. However, distinguishing genes containing driving mutations rather than simply passenger SSVs from a cohort sequenced cancer samples requires sophisticated computational approaches.
20/20+ integrates many features indicative of positive selection to predict oncogenes and tumor suppressor genes from small somatic variants. 
The features capture mutational clustering, conservation, mutation *in silico* pathogenicity scores, mutation consequence types, protein interaction network connectivity, and other covariates (e.g. replication timing).
Contrary to methods based on mutation rate, 20/20+ uses ratiometric features of mutations by normalizing for the total number of mutations in a gene. This decouples the genes from gene-level differences in background mutation rate. 20/20+ uses monte carlo simulations to evaluate the significance of random forest scores based on an estimated p-value from an empirical null distribution.

## Documentation

[![Documentation Status](http://readthedocs.org/projects/2020plus/badge/?version=latest)](http://2020plus.readthedocs.io/en/latest/?badge=latest)

Please see the [documentation](http://2020plus.readthedocs.io/) on readthedocs.

## Releases

You can download [releases](https://github.com/KarchinLab/2020plus/releases) on github.

## Installation

[![Build Status](https://travis-ci.org/KarchinLab/2020plus.svg?branch=master)](https://travis-ci.org/KarchinLab/2020plus)

We recommend that you install the dependencies for 20/20+ through [conda](https://conda.io/miniconda.html). Once conda is installed, setting up the environment is done as follows:

```bash
$ conda env create -f environment_python.yml  # install dependencies for python
$ source activate 2020plus  # activate the 20/20+ conda environment
$ conda install r r-randomForest rpy2  # install the R related dependencies
```

Every time you wish to run 20/20+, you will then need to activate the "2020plus" conda environment.

```bash
$ source activate 2020plus
```

The 20/20+ conda environment can also be deactivated.

```bash
$ source deactivate 2020plus
```
