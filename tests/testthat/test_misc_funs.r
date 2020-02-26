### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

source("tolerances.r")

# unsymmetric matrix
unsym_mat <- matrix(0.5, 2, 2); diag(unsym_mat) <- 1; unsym_mat[1, 2] <- 0.2

# matrix with a missing value
mat_w_mis <- matrix(0.5, 2, 2); diag(mat_w_mis) <- 1; mat_w_mis[1, 1] <- NA

# negative-definite matrix
neg_def_mat <- matrix(-0.8, 3, 3); diag(neg_def_mat) <- 1

# matrix as a data frame
dat_fra_mat <- matrix(0.8, 3, 3); diag(dat_fra_mat) <- 1
dat_fra_mat <- as.data.frame(dat_fra_mat)

# matrix with values out of bounds
mat_out_bou <- matrix(1.5, 2, 2); diag(unsym_mat) <- 1

# matrix with diagonal values other than 1
mat_out_diag <- matrix(0.9, 2, 2)

# an appropriate matrix (to test the dimensions with the vector of p-values)
approp_mat <- matrix(0.5, 2, 2); diag(approp_mat) <- 1 

context("Checking errors")

test_that("Errors are thrown correctly.", {
  
  expect_error(fisher(), "Argument 'p' must be specified.")
  expect_error(fisher(matrix("a", 2, 2)), "Argument 'p' must be a numeric vector.")
  expect_error(fisher(c(0.1, NA)), "Values in 'p' vector must not contain NAs.")
  expect_error(fisher(c(1.1, 0.1)), "Values in 'p' vector \\(i.e., the p-values\\) must be between 0 and 1.")
  
  expect_error(mvnconv(unsym_mat, target = "m2lp"), "Argument 'R' must be a symmetric matrix.")
  expect_error(meff(mat_w_mis, method = "liji"), "Values in 'R' vector must not contain NAs.")
  expect_warning(empirical(neg_def_mat, method = "fisher"), "Matrix 'R' is not positive definite. Used Matrix::nearPD\\(\\) to make 'R' positive definite.")
  expect_error(meff(mat_out_bou, method = "liji"), "Argument 'R' must be a correlation matrix, but contains values outside \\[-1,1\\].")
  
  expect_error(fisher(runif(3), adjust = "liji", R = approp_mat), "Length of 'p' vector \\(3\\) does not match the dimensions of the 'R' matrix \\(2,2\\).")
  expect_warning(fisher(runif(2), R = approp_mat))
  expect_warning(fisher(runif(2), m = 3, R = approp_mat))
  expect_warning(fisher(runif(2), m = 3))
  
  expect_error(empirical(approp_mat, method = "fisher", side = c(1, 2)), "Argument 'side' must be of length 1.")
  expect_error(empirical(approp_mat, method = "fisher", side = 3), "Argument 'side' must be either 1 or 2.")
  
  expect_warning(fisher(runif(2), adjust = "empirical", R = approp_mat, emp.dist = runif(10), threshold = 0.5))

  expect_error(fisher(runif(2), adjust = "empirical", R = approp_mat, size = "a"), "Argument 'size' must be numeric. See help\\(fisher\\).")
  expect_error(fisher(runif(2), adjust = "empirical", R = approp_mat, size = -1), "Values in 'size' must be >= 1. See help\\(fisher\\).")
  
  expect_warning(fisher(runif(2), adjust = "empirical", R = approp_mat, threshold = 0.5))
  
  expect_error(fisher(runif(2), adjust = "empirical", R = approp_mat, size = c(100, 1000)), "Argument 'threshold' must be specified when 'size' is a vector. See help\\(fisher\\).")
  expect_error(fisher(runif(2), adjust = "empirical", R = approp_mat, size = c(100, 1000), threshold = "a"), "Argument 'threshold' must be numeric. See help\\(fisher\\).")
  expect_error(fisher(runif(2), adjust = "empirical", R = approp_mat, size = c(100, 1000), threshold = 1.1), "Values in 'threshold' must be between 0 and 1. See help\\(fisher\\).")
  expect_error(fisher(runif(2), adjust = "empirical", R = approp_mat, size = c(100, 1000), threshold = c(0.3, 0.3, 0.1)), "Length of 'threshold' argument is not compatible with length of 'size' argument. See help\\(fisher\\).")
  expect_error(fisher(runif(2), adjust = "empirical", R = approp_mat, size = c(100, 1000), threshold = c(0.3, 0.3, 0.1)), "Length of 'threshold' argument is not compatible with length of 'size' argument. See help\\(fisher\\).")
  
  out <- capture.output(fisher(runif(2), adjust = "empirical", R = approp_mat, size = c(100, 1000), threshold = c(0.3), verbose = TRUE))
  
  expect_error(fisher(runif(2), adjust = "empirical"), "Argument 'R' must be specified when using an adjustment method.")
  expect_error(stouffer(runif(2), adjust = "empirical"), "Argument 'R' must be specified when using an adjustment method.")
  expect_error(invchisq(runif(2), adjust = "empirical"), "Argument 'R' must be specified when using an adjustment method.")
  expect_error(binotest(runif(2), adjust = "empirical"), "Argument 'R' must be specified when using an adjustment method.")
  expect_error(bonferroni(runif(2), adjust = "empirical"), "Argument 'R' must be specified when using an adjustment method.")
  expect_error(tippett(runif(2), adjust = "empirical"), "Argument 'R' must be specified when using an adjustment method.")
  
  expect_error(fisher(runif(2), adjust = "liji", R = mat_out_diag), "Diagonal values in 'R' must all be equal to 1.")
  
  expect_warning(fisher(runif(2), adjust = "empirical", R = approp_mat, size = c(100, 1000, 10000), threshold = rep(0.5, 3)))
  
  expect_warning(fisher(runif(2), adjust = "generalized", R = mvnconv(approp_mat, target = "z")))
  
})

test_that("Conversions work correctly.", {
  
  meff_neg_def_mat <- meff(neg_def_mat, method = "liji")
  meff_dat_fra_mat <- meff(dat_fra_mat, method = "liji")
  
  expect_equivalent(meff_neg_def_mat, 4, tolerance = m_tol)
  expect_equivalent(meff_dat_fra_mat, 2, tolerance = m_tol)
  
  set.seed(1234)
  meff_nearpd <- fisher(runif(3), adjust = "liji", R = nearPD(neg_def_mat, corr = TRUE)$mat)
  expect_equivalent(c(meff_nearpd$p), 0.3917173, tolerance = p_tol)

})
