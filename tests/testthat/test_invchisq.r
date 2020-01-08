### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

source("tolerances.r")

context("Checking invchisq() function")

test_that("invchisq() works correctly under independence.", {

  res <- invchisq(grid2ip.p)

  expect_equivalent(c(res$p), 6.915992e-08, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 77.9029090, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 23, tolerance = df_tol)

})

test_that("invchisq() works correctly with effective number of tests.", {

  res_nyh <- invchisq(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj  <- invchisq(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- invchisq(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- invchisq(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))

  expect_equivalent(c(res_nyh$p), 0.0000001, tolerance = p_tol)
  expect_equivalent(c(res_nyh$statistic), 74.5158260, tolerance = stat_tol)

  expect_equivalent(c(res_lj$p), 0.0000002, tolerance = p_tol)
  expect_equivalent(c(res_lj$statistic), 71.1287430, tolerance = stat_tol)

  expect_equivalent(c(res_gao$p), 6.915992e-08, tolerance = p_tol)
  expect_equivalent(c(res_gao$statistic), 77.9029090, tolerance = stat_tol)

  expect_equivalent(c(res_gal$p), 0.0000004, tolerance = p_tol)
  expect_equivalent(c(res_gal$statistic), 67.7416600, tolerance = stat_tol)

})

test_that("invchisq() works correctly with empirically-derived null distributions.", {

  set.seed(1234)
  res <- invchisq(grid2ip.p, adjust = "empirical", R = grid2ip.ld)

  expect_equivalent(c(res$p), 0.0024998, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 77.9029090, tolerance = stat_tol)
  expect_equivalent(c(res$ci[1]), 0.0016183, tolerance = stat_tol)
  expect_equivalent(c(res$ci[2]), 0.0036879, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 23, tolerance = df_tol)

  set.seed(1234)
  res <- invchisq(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)

  expect_equivalent(c(res$p), 0.0022900, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 77.9029090, tolerance = stat_tol)
  expect_equivalent(c(res$ci[1]), 0.0020032, tolerance = stat_tol)
  expect_equivalent(c(res$ci[2]), 0.0026062, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 23, tolerance = df_tol)

  set.seed(1234)
  res <- invchisq(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)

  expect_equivalent(c(res$p), 0.0023240, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 77.9029090, tolerance = stat_tol)
  expect_equivalent(c(res$ci[1]), 0.0022306, tolerance = stat_tol)
  expect_equivalent(c(res$ci[2]), 0.0024203, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 23, tolerance = df_tol)

  set.seed(1234)
  res <- invchisq(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))

  expect_equivalent(c(res$p), 0.0026900, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 77.9029090, tolerance = stat_tol)
  expect_equivalent(c(res$ci[1]), 0.0023785, tolerance = stat_tol)
  expect_equivalent(c(res$ci[2]), 0.0030309, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 23, tolerance = df_tol)

})

test_that("invchisq() works correctly under multivariate theory.", {

  res1 <- invchisq(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 1))

  expect_equivalent(c(res1$p), 0.0000845, tolerance = p_tol)
  expect_equivalent(c(res1$statistic), 38.4072406, tolerance = stat_tol)
  expect_equivalent(attributes(res1$statistic)$df, 11.3393267, tolerance = df_tol)

  res2 <- invchisq(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 2))

  expect_equivalent(c(res2$p), 0.0010150, tolerance = p_tol)
  expect_equivalent(c(res2$statistic), 24.9527814, tolerance = stat_tol)
  expect_equivalent(attributes(res2$statistic)$df, 7.3670416, tolerance = df_tol)

})
