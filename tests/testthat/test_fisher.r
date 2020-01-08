### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

source("tolerances.r")

context("Checking fisher() function")

test_that("fisher() works correctly under independence.", {

  res <- fisher(grid2ip.p)

  expect_equivalent(c(res$p), 2.757113e-08, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 118.2920513, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

})

test_that("fisher() works correctly with effective number of tests.", {

  res_nyh <- fisher(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj  <- fisher(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- fisher(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- fisher(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))

  expect_equivalent(c(res_nyh$p), 5.267347e-08, tolerance = p_tol)
  expect_equivalent(c(res_nyh$statistic), 113.1489187, tolerance = stat_tol)

  expect_equivalent(c(res_lj$p), 0.0000001, tolerance = p_tol)
  expect_equivalent(c(res_lj$statistic), 108.0057860, tolerance = stat_tol)

  expect_equivalent(c(res_gao$p), 2.757113e-08, tolerance = p_tol)
  expect_equivalent(c(res_gao$statistic), 118.2920513, tolerance = stat_tol)

  expect_equivalent(c(res_gal$p), 0.0000002, tolerance = p_tol)
  expect_equivalent(c(res_gal$statistic), 102.8626533, tolerance = stat_tol)

})

test_that("fisher() works correctly with empirically-derived null distributions.", {

  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld)

  expect_equivalent(c(res$p), 0.0020998, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 118.2920513, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0013003, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0032080, tolerance = stat_tol * emp_sca)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)

  expect_equivalent(c(res$p), 0.0020300, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 118.2920513, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0017605, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0023289, tolerance = stat_tol * emp_sca)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)

  expect_equivalent(c(res$p), 0.0020410, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 118.2920513, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0019535, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0021314, tolerance = stat_tol * emp_sca)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))

  expect_equivalent(c(res$p), 0.0024000, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 118.2920513, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0021062, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0027232, tolerance = stat_tol * emp_sca)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

})

test_that("fisher() works correctly under multivariate theory.", {

  res1 <- fisher(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 1))

  expect_equivalent(c(res1$p), 0.0000312, tolerance = p_tol)
  expect_equivalent(c(res1$statistic), 62.9090102, tolerance = stat_tol)
  expect_equivalent(attributes(res1$statistic)$df, 24.4633045, tolerance = df_tol)

  res2 <- fisher(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 2))

  expect_equivalent(c(res2$p), 0.0007643, tolerance = p_tol)
  expect_equivalent(c(res2$statistic), 38.3482083, tolerance = stat_tol)
  expect_equivalent(attributes(res2$statistic)$df, 14.9123932, tolerance = df_tol)

})
