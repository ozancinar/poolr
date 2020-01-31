### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

source("tolerances.r")

context("Checking mvnconv() function")

test_that("empirical() works correctly.", {
  
  set.seed(1234)
  emp_test_alpha <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, alpha = 0.1)
  
  set.seed(1234)
  emp_test_batch <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000, batchsize = 300)
  
  set.seed(1234)
  emp_test_side1 <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000, side = 1)
  
  expect_equivalent(emp_test_alpha$p, 0.00079992, tolerance = p_tol)
  expect_equivalent(emp_test_batch$p, 0.002997003, tolerance = p_tol)
  expect_equivalent(emp_test_side1$p, 0.000999001, tolerance = p_tol)
  
})

test_that("The arguments of empirical() are checked correctly.", {
  
  expect_error(empirical(method = "fisher"), "Argument 'R' must be specified.")
  expect_error(fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100, batchsize = 1000))
  # expect_error(fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000000, side = 1), "Matrix to be generated is too large. Try setting 'batchsize' \\(or to a lower number if it was set\\).")
  
})
