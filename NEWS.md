# poolr 0.9-9 (2020-06-11)

- added a warning message to `meff(method = "galwey")` if there are negative eigenvalues

- `mvnconv()` now checks if `R` is non-negative definite

# poolr 0.9-8 (2020-05-08)

- the HTML help files now show rendered equations with the help of the `mathjaxr` package

- implemented a simplified version of `Matrix::nearPD()`

- dependence on the package `Matrix` was removed

# poolr 0.9-6 (2020-03-18)

- added a more specific test on `p` and `eigen` that they are `numeric` vectors

# poolr 0.9-3 (2020-02-26)

- improved the `pkgdown` docs and added a quick start guide

- changed the way the pseudo replicates are generated in `empirical()` to a more stable method

- slight improvements to the output of `print.poolr()` when using the effective number of tests or empirical distribution adjustments

- `mvnconv()` now uses the variances from the lookup table instead of `cov2cor()` for the transformation when `cov2cor=TRUE`

- added a check on `R` (where appropriate) that its diagonal values are all equal to 1

- added a check on `p` to convert it into a `numeric` vector if it is a `matrix` with 1 row

# poolr 0.8-2 (2020-02-12)

- first version for CRAN
