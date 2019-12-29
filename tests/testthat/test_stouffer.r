
### library(library); library(testthat); Sys.setenv(NOT_CRAN="true")

p_tol <- 1e-06  # user-defined tolerance (for significant p-values)
stat_tol <- 1e-04  # user-defined tolerance (for test statistics and ci bounds)

context("Checking stouffer() function")

test_that("stouffer() works correctly under independence.", {
  
  res <- stouffer(grid2ip.p)

  expect_equivalent(as.numeric(res$p), 3.62e-08, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 5.3851, tolerance = stat_tol)
  
})


test_that("stouffer() works correctly with effective number of tests.", {
  
  res_nyh <- stouffer(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj <- stouffer(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- stouffer(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- stouffer(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  
  expect_equivalent(as.numeric(res_nyh$p), 6.94e-08, tolerance = p_tol)
  expect_equivalent(as.numeric(res_nyh$statistic), 5.2668, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_lj$p), 1.33e-07, tolerance = p_tol)
  expect_equivalent(as.numeric(res_lj$statistic), 5.1457, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_gao$p), 3.62e-08, tolerance = p_tol)
  expect_equivalent(as.numeric(res_gao$statistic), 5.3851, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_gal$p), 2.56e-07, tolerance = p_tol)
  expect_equivalent(as.numeric(res_gal$statistic), 5.0217, tolerance = stat_tol)
  
})


test_that("stouffer() works correctly with empirically-derived null distributions.", {
  
  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld)
  
  expect_equivalent(as.numeric(res$p), 0.00130, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 5.3851, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0007, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0022, tolerance = stat_tol)
  
  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)
  
  expect_equivalent(as.numeric(res$p), 0.00159, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 5.3851, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0014, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0019, tolerance = stat_tol)
  
  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)
  
  expect_equivalent(as.numeric(res$p), 0.001672, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 5.3851, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0016, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0018, tolerance = stat_tol)
  
  set.seed(1234)
  res <- stouffer(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))
  
  expect_equivalent(as.numeric(res$p), 0.00194, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 5.3851, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0017, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0022, tolerance = stat_tol)
  
})


test_that("stouffer() works correctly under multivariate theory.", {
  
  res1 <- stouffer(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 1))
  
  expect_equivalent(as.numeric(res1$p), 1.13e-05, tolerance = p_tol)
  expect_equivalent(as.numeric(res1$statistic), 4.2366, tolerance = stat_tol)
  
  res2 <- stouffer(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 2))
  
  expect_equivalent(as.numeric(res2$p), 0.000402, tolerance = p_tol)
  expect_equivalent(as.numeric(res2$statistic), 3.3517, tolerance = stat_tol)
  
})
