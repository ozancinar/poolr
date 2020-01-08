poolr: Package for Pooling the Results from (Dependent) Tests
=============================================================

[![Build Status](https://travis-ci.org/ozancinar/poolr.svg?branch=master)](https://travis-ci.org/ozancinar/poolr)
![CRAN Version](https://www.r-pkg.org/badges/version/poolr)
![devel Version](https://img.shields.io/badge/devel-0.1--8-brightgreen.svg)

## Description

The `poolr` package contains functions for pooling/combining the results (i.e., p-values) from (dependent) hypothesis tests. Included are Fisher's method, Stouffer's method, the inverse chi-square method, the Bonferroni method, Tippett's method, and the binomial test. Each method can be adjusted based on an estimate of the effective number of tests or using empirically derived null distribution using pseudo replicates. For Fisher's, Stouffer's, and the inverse chi-square method, direct generalizations based on multivariate theory are also available (leading to Brown's method, Strube's method, and the generalized inverse chi-square method).

## Documentation

You can read the documentation of the `poolr` package online at [https://ozancinar.github.io/poolr/](https://ozancinar.github.io/poolr/) (where it is nicely formatted, equations are shown correctly, and the output from all examples is provided).

## Installation

After installing the [remotes](https://cran.r-project.org/package=remotes) package with ```install.packages("remotes")```, the development version of the `poolr` package can be installed with:
```r
remotes::install_github("ozancinar/poolr")
```

## Meta

The poolr package was written by Ozan Cinar and [Wolfgang Viechtbauer](http://www.wvbauer.com/). It is licensed under the [GNU General Public License Version 2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.txt). To report any issues or bugs, please go [here](https://github.com/ozancinar/poolr/issues).
