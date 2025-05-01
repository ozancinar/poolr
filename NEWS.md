# poolr 1.2-0 (2025-05-01)

- added the reference to the publication in the Journal of Statistical Software

- fixed an issue where the Bonferroni and Tippett methods were not identical with `adjust = "empirical"` when the observed Bonferroni-adjusted p-value is equal to 1

- added the method by Chen & Liu (2011) for estimating the effective number of tests via `adjust = "chen"`

# poolr 1.0-0 (2021-06-02)

- changed name of `binotest()` function to `binomtest()`

- the HTML help files now show rendered equations with the help of the `mathjaxr` package

- increased resolution of `mvnlookup` table (now in steps of .001)

- `meff()` function now issues a warning if there are negative eigenvalues (and if they were set to 0 for `method="galway"`)

- added `nearpd` argument to all base functions; if `TRUE`, a negative definite `R` matrix will be turned into the nearest positive semi-definite matrix (only for `adjust="empirical"` and `adjust="generalized"`)

- implemented a simplified version of `Matrix::nearPD()`; hence, dependence on the package `Matrix` was removed

- added a more specific test on `p` and `eigen` that they are `numeric` vectors

- improved the `pkgdown` docs and added a quick start guide

- changed the way the pseudo replicates are generated in `empirical()` to a more stable method

- slight improvements to the output of `print.poolr()` when using the effective number of tests or empirical distribution adjustments

- `mvnconv()` now uses the variances from the lookup table instead of `cov2cor()` for the transformation when `cov2cor=TRUE`

- added a check on `R` (where appropriate) that its diagonal values are all equal to 1

- added a check on `p` to convert it into a `numeric` vector if it is a `matrix` with 1 row

# poolr 0.8-2 (2020-02-12)

- first version for CRAN
