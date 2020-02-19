### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

source("tolerances.r")

context("Checking mvnconv() function")

test_that("mvnconv() works correctly.", {
  
  mvnconv_vec_ex1 <- mvnconv(c(0.1, 0.2, 0.3), target = "m2lp")
  mvnconv_vec_ex2 <- mvnconv(c(0.1, 0.2, 0.3), target = "m2lp", cov2cor = TRUE)

  expect_equivalent(mvnconv_vec_ex1, c(0.0390, 0.1563, 0.3519), tolerance = stat_tol)
  expect_equivalent(mvnconv_vec_ex2, c(0.00975, 0.039075, 0.087975), tolerance = stat_tol)
  
})

test_that("The arguments of mvnconv() are checked correctly.", {
  
  expect_error(mvnconv(target = "m2lp"), "Argument 'R' must be specified.")
  expect_error(mvnconv(grid2ip.ld), "Argument 'target' must be specified.")
  
})
