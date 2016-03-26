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
* pandas
* scikit-learn
* rpy2

To install these packages via `pip` you can use the following command:

```bash
$ pip install -r requirements.txt
```
