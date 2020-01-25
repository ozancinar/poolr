### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

source("tolerances.r")

context("Checking binotest() function")

test_that("binotest() works correctly under independence.", {

  res <- binotest(grid2ip.p)

  expect_equivalent(c(res$p), 3.763872e-09, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 11, tolerance = stat_tol)

})

test_that("binotest() works correctly with effective number of tests.", {

  res_nyh <- binotest(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj  <- binotest(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- binotest(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- binotest(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))

  expect_equivalent(c(res_nyh$p), 2.057712e-09, tolerance = p_tol)
  expect_equivalent(c(res_nyh$statistic), 11, tolerance = stat_tol)

  expect_equivalent(c(res_lj$p), 2.067037e-08, tolerance = p_tol)
  expect_equivalent(c(res_lj$statistic), 11, tolerance = stat_tol)

  expect_equivalent(c(res_gao$p), 3.763872e-09, tolerance = p_tol)
  expect_equivalent(c(res_gao$statistic), 11, tolerance = stat_tol)

  expect_equivalent(c(res_gal$p), 1.134072e-08, tolerance = p_tol)
  expect_equivalent(c(res_gal$statistic), 11, tolerance = stat_tol)

})

test_that("binotest() works correctly with empirically-derived null distributions.", {

  set.seed(1234)
  res <- binotest(grid2ip.p, adjust = "empirical", R = grid2ip.ld)

  expect_equivalent(c(res$p), 0.00059994, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 11, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0002201982, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.001305356, tolerance = stat_tol * emp_sca)

  set.seed(1234)
  res <- binotest(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)

  expect_equivalent(c(res$p), 0.0003899961, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 11, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0002773396, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0005330993, tolerance = stat_tol * emp_sca)

  set.seed(1234)
  res <- binotest(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)

  expect_equivalent(c(res$p), 0.0004289996, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 11, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.000389367, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0004715709, tolerance = stat_tol * emp_sca)

  set.seed(1234)
  res <- binotest(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))

  expect_equivalent(c(res$p), 0.0004799952, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 11, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0003539311, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0006363542, tolerance = stat_tol * emp_sca)

})
