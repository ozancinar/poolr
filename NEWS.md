# poolr 0.9-1 (2020-02-20)

- changed the way the pseudo replicates are generated in `empirical()` to a more stable method

- slight improvements to the output of `print.poolr()` when using the effective number of tests or empirical distribution adjustments

- `mvnconv()` now uses the variances from the lookup table instead of `cov2cor()` for the transformation when `cov2cor=TRUE`

- added a check on `R` (where appropriate) that its diagonal values are all equal to 1

# poolr 0.8-2 (2020-02-12)

- first version for CRAN
