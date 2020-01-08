### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

source("tolerances.r")

context("Checking stouffer() function")

test_that("stouffer() works correctly under independence.", {

  res <- stouffer(grid2ip.p)

  expect_equivalent(c(res$p), 3.619282e-08, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 5.3851466, tolerance = stat_tol)

})

test_that("stouffer() works correctly with effective number of tests.", {

  res_nyh <- stouffer(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj  <- stouffer(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- stouffer(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- stouffer(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))

  expect_equivalent(c(res_nyh$p), 6.941975e-08, tolerance = p_tol)
  expect_equivalent(c(res_nyh$statistic), 5.2667773, tolerance = stat_tol)

  expect_equivalent(c(res_lj$p), 0.0000001, tolerance = p_tol)
  expect_equivalent(c(res_lj$statistic), 5.1456857, tolerance = stat_tol)

  expect_equivalent(c(res_gao$p), 3.619282e-08, tolerance = p_tol)
  expect_equivalent(c(res_gao$statistic), 5.3851466, tolerance = stat_tol)

  expect_equivalent(c(res_gal$p), 0.0000003, tolerance = p_tol)
  expect_equivalent(c(res_gal$statistic), 5.0216751, tolerance = stat_tol)

})

test_that("stouffer() works correctly with empirically-derived null distributions.", {

  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld)

  expect_equivalent(c(res$p), 0.0012999, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 5.3851466, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0006923, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0022218, tolerance = stat_tol * emp_sca)

  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)

  expect_equivalent(c(res$p), 0.0015900, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 5.3851466, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0013526, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0018570, tolerance = stat_tol * emp_sca)

  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)

  expect_equivalent(c(res$p), 0.0016720, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 5.3851466, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0015929, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0017540, tolerance = stat_tol * emp_sca)

  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))

  expect_equivalent(c(res$p), 0.0019400, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 5.3851466, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0016768, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.0022327, tolerance = stat_tol * emp_sca)

})

test_that("stouffer() works correctly under multivariate theory.", {

  res1 <- stouffer(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 1))

  expect_equivalent(c(res1$p), 0.0000113, tolerance = p_tol)
  expect_equivalent(c(res1$statistic), 4.2366598, tolerance = stat_tol)

  res2 <- stouffer(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 2))

  expect_equivalent(c(res2$p), 0.0004015, tolerance = p_tol)
  expect_equivalent(c(res2$statistic), 3.3517422, tolerance = stat_tol)

})
