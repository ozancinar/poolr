### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

source("tolerances.r")

context("Checking invchisq() function")

test_that("invchisq() works correctly under independence.", {

  res <- invchisq(grid2ip.p)

  expect_equivalent(c(res$p), 6.92e-08, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 77.9029, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 23, tolerance = df_tol)

})

test_that("invchisq() works correctly with effective number of tests.", {

  res_nyh <- invchisq(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj  <- invchisq(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- invchisq(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- invchisq(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))

  expect_equivalent(c(res_nyh$p), 1.26e-07, tolerance = p_tol)
  expect_equivalent(c(res_nyh$statistic), 74.5158, tolerance = stat_tol)

  expect_equivalent(c(res_lj$p), 2.31e-07, tolerance = p_tol)
  expect_equivalent(c(res_lj$statistic), 71.1287, tolerance = stat_tol)

  expect_equivalent(c(res_gao$p), 6.92e-08, tolerance = p_tol)
  expect_equivalent(c(res_gao$statistic), 77.9029, tolerance = stat_tol)

  expect_equivalent(c(res_gal$p), 4.24e-07, tolerance = p_tol)
  expect_equivalent(c(res_gal$statistic), 67.7417, tolerance = stat_tol)

})

test_that("invchisq() works correctly with empirically-derived null distributions.", {

  set.seed(1234)
  res <- invchisq(grid2ip.p, adjust = "empirical", R = grid2ip.ld)

  expect_equivalent(c(res$p), 0.00250, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 77.9029, tolerance = stat_tol)
  expect_equivalent(c(res$ci[1]), 0.0016, tolerance = stat_tol)
  expect_equivalent(c(res$ci[2]), 0.0037, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 23, tolerance = df_tol)

  set.seed(1234)
  res <- invchisq(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)

  expect_equivalent(c(res$p), 0.00229, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 77.9029, tolerance = stat_tol)
  expect_equivalent(c(res$ci[1]), 0.0020, tolerance = stat_tol)
  expect_equivalent(c(res$ci[2]), 0.0026, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 23, tolerance = df_tol)

  set.seed(1234)
  res <- invchisq(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)

  expect_equivalent(c(res$p), 0.002324, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 77.9029, tolerance = stat_tol)
  expect_equivalent(c(res$ci[1]), 0.0022, tolerance = stat_tol)
  expect_equivalent(c(res$ci[2]), 0.0024, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 23, tolerance = df_tol)

  set.seed(1234)
  res <- invchisq(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))

  expect_equivalent(c(res$p), 0.00269, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 77.9029, tolerance = stat_tol)
  expect_equivalent(c(res$ci[1]), 0.0024, tolerance = stat_tol)
  expect_equivalent(c(res$ci[2]), 0.0030, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 23, tolerance = df_tol)

})

test_that("invchisq() works correctly under multivariate theory.", {

  res1 <- invchisq(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 1))

  expect_equivalent(c(res1$p), 8.45e-05, tolerance = p_tol)
  expect_equivalent(c(res1$statistic), 38.4072, tolerance = stat_tol)
  expect_equivalent(attributes(res1$statistic)$df, 11.34, tolerance = df_tol)

  res2 <- invchisq(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 2))

  expect_equivalent(c(res2$p), 0.001015, tolerance = p_tol)
  expect_equivalent(c(res2$statistic), 24.9528, tolerance = stat_tol)
  expect_equivalent(attributes(res2$statistic)$df, 7.37, tolerance = df_tol)

})
