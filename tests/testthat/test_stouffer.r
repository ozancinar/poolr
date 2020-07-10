### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")

source("tolerances.r")

context("Checking stouffer() function")

test_that("stouffer() works correctly under independence.", {

  res <- stouffer(grid2ip.p)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 1.655433e-09, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 5.915392, tolerance = stat_tol)

})

test_that("stouffer() works correctly with effective number of tests.", {

  res_nyh <- stouffer(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj  <- stouffer(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- stouffer(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- stouffer(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_user <- stouffer(grid2ip.p, m = 18)

  out <- capture.output(print(res_nyh))
  out <- capture.output(print(res_lj))
  out <- capture.output(print(res_gao))
  out <- capture.output(print(res_gal))
  out <- capture.output(print(res_user))

  expect_equivalent(c(res_nyh$p), 3.617704e-09, tolerance = p_tol)
  expect_equivalent(c(res_nyh$statistic), 5.785367, tolerance = stat_tol)

  expect_equivalent(c(res_lj$p), 7.913328e-09, tolerance = p_tol)
  expect_equivalent(c(res_lj$statistic), 5.652353, tolerance = stat_tol)

  expect_equivalent(c(res_gao$p), 1.655433e-09, tolerance = p_tol)
  expect_equivalent(c(res_gao$statistic), 5.915392, tolerance = stat_tol)

  expect_equivalent(c(res_gal$p), 1.732717e-08, tolerance = p_tol)
  expect_equivalent(c(res_gal$statistic), 5.516131, tolerance = stat_tol)

  expect_equivalent(c(res_user$p), 8.336258e-08, tolerance = p_tol)
  expect_equivalent(c(res_user$statistic), 5.233062, tolerance = stat_tol)

})

test_that("stouffer() works correctly with empirically-derived null distributions.", {

  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.00049995, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 5.915392, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0001623517, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.001166328, tolerance = p_tol * emp_sca)

  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.0007599924, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 5.915392, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0005988329, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0009511528, tolerance = p_tol * emp_sca)

  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.0007119993, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 5.915392, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0006625978, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0007683274, tolerance = p_tol * emp_sca)

  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.0007899921, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 5.915392, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0006254925, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0009844701, tolerance = p_tol * emp_sca)

})

test_that("stouffer() works correctly under multivariate theory.", {

  res1 <- stouffer(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 1))
  out <- capture.output(print(res1))

  expect_equivalent(c(res1$p), 1.67123e-06, tolerance = p_tol)
  expect_equivalent(c(res1$statistic), 4.648569, tolerance = stat_tol)

  res2 <- stouffer(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 2))
  out <- capture.output(print(res2))

  expect_equivalent(c(res2$p), 0.000112044, tolerance = p_tol)
  expect_equivalent(c(res2$statistic), 3.690188, tolerance = stat_tol)

})
