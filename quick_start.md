# Quick Start Guide

The package contains six 'base methods' for pooling p-values:

* `fisher()`: for Fisher's method,
* `stouffer()`: for Stouffer's method,
* `invchisq()`: for the inverse chi-square method,
* `binomtest()`: for the binomial test,
* `bonferroni()`: for the Bonferroni method,
* `tippett()`: for Tippett's method.

For example, we can combine a set of independent p-values using [Fisher's method](https://en.wikipedia.org/wiki/Fisher's_method) as follows:

```r
library(poolr)
pvals <- c(0.02, 0.03, 0.08, 0.20)
fisher(pvals)
```
```
combined p-values with:      Fisher's method
number of p-values combined: 4
test statistic:              23.107 ~ chi-square(8)
adjustment:                  none
combined p-value:            0.003228942
```

More interesting are cases where the p-values are not independent. For example,

```r
round(grid2ip.p[1:5], digits = 5)
```
```
 rs10267908 rs112305062 rs117541653  rs11761490  rs11773436
    0.01137     0.50636     0.12303     0.09992     0.00169
```

shows the first 5 p-values from (two-sided) tests of the association between 23 [single-nucleotide polymorphism](https://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) (SNP) in the *GRID2IP* gene and depressive symptoms. Due to [linkage disequilibrium](https://en.wikipedia.org/wiki/Linkage_disequilibrium) (LD), the SNPs are not independent and hence neither will the p-values. The following shows the first 5 rows and columns of the LD correlation matrix:

```r
round(grid2ip.ld[1:5,1:5], digits = 3)
```
```
            rs10267908 rs112305062 rs117541653 rs11761490 rs11773436
rs10267908       1.000       0.187      -0.192     -0.130     -0.389
rs112305062      0.187       1.000       0.145     -0.010     -0.260
rs117541653     -0.192       0.145       1.000     -0.097      0.099
rs11761490      -0.130      -0.010      -0.097      1.000     -0.019
rs11773436      -0.389      -0.260       0.099     -0.019      1.000
```

We can adjust the various base methods to account for the dependence using this correlation matrix. For example, we can use an estimate of the effective number of tests based on Li and Ji (2005) to adjust the test statistic of Fisher's method as a way to account for the dependence:

```{r}
fisher(grid2ip.p, adjust = "liji", R = grid2ip.ld)
```
```
combined p-values with:      Fisher's method
number of p-values combined: 23
test statistic:              83.14 ~ chi-square(30)
adjustment:                  effective number of tests (m = 15; Li & Ji, 2005)
combined p-value:            6.92e-07
```

Alternatively, we can use the generalization of Fisher's method described by Brown (1975) to combine the p-values:

```{r}
fisher(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld))
```
```
combined p-values with:      Fisher's method
number of p-values combined: 23
test statistic:              41.532 ~ chi-square(14.986)
adjustment:                  Brown's method
combined p-value:            0.000263
```

Finally, one can empirically obtain the null distribution of Fisher's method using pseudo replicates and compute the combined p-value based on that (which closely approximates a 'proper' permutation test, but runs in a fraction of the time):

```{r}
set.seed(123)
fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld)
```
```
combined p-values with:      Fisher's method
number of p-values combined: 23
test statistic:              127.482 ~ chi-square(46)
adjustment:                  empirical distribution (size = 10000)
combined p-value:            0.0012 (95% CI: 0.00062, 0.0021)
```

Since this method is stochastic, we manually specify the seed for the random number generator for reproducibility.

The examples above cover only part of the functionality of the package. You can read the documentation of all functions [here](https://ozancinar.github.io/poolr/reference/index.html).
