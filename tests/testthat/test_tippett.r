### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

source("tolerances.r")

context("Checking tippett() function")

test_that("tippett() works correctly under independence.", {

  res <- tippett(grid2ip.p)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.03810371, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 0.001687646, tolerance = stat_tol)

})

test_that("tippett() works correctly with effective number of tests.", {

  res_nyh <- tippett(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj  <- tippett(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- tippett(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- tippett(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  
  out <- capture.output(print(res_nyh))
  out <- capture.output(print(res_lj))
  out <- capture.output(print(res_gao))
  out <- capture.output(print(res_gal))

  expect_equivalent(c(res_nyh$p), 0.03647763, tolerance = p_tol)
  expect_equivalent(c(res_nyh$statistic), 0.001687646, tolerance = stat_tol)

  expect_equivalent(c(res_lj$p), 0.03484879, tolerance = p_tol)
  expect_equivalent(c(res_lj$statistic), 0.001687646, tolerance = stat_tol)

  expect_equivalent(c(res_gao$p), 0.03810371, tolerance = p_tol)
  expect_equivalent(c(res_gao$statistic), 0.001687646, tolerance = stat_tol)

  expect_equivalent(c(res_gal$p), 0.03321721, tolerance = p_tol)
  expect_equivalent(c(res_gal$statistic), 0.001687646, tolerance = stat_tol)

})

test_that("tippett() works correctly with empirically-derived null distributions.", {

  set.seed(1234)
  res <- tippett(grid2ip.p, adjust = "empirical", R = grid2ip.ld)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.03059694, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 0.001687646, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.02730888, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.03416195, tolerance = p_tol * emp_sca)

  set.seed(1234)
  res <- tippett(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.03067969, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 0.001687646, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0296195, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.03176722, tolerance = p_tol * emp_sca)

  set.seed(1234)
  res <- tippett(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.03017697, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 0.001687646, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.02984253, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.03051413, tolerance = p_tol * emp_sca)

  set.seed(1234)
  res <- tippett(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.0309969, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 0.001687646, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.02768743, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.03458305, tolerance = p_tol * emp_sca)

})
