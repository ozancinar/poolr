
### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

# p_tol <- 1e-06  # user-defined tolerance (for significant p-values)
# stat_tol <- 1e-04  # user-defined tolerance (for test statistics and ci bounds)
# df_tol <- 1e-02  # user-defined tolerance (for df)

p_tol <- 0.1
stat_tol <- 0.1
df_tol <- 0.1

context("Checking fisher() function")

test_that("fisher() works correctly under independence.", {
  
  res <- fisher(grid2ip.p)

  expect_equivalent(as.numeric(res$p), 2.75e-08, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 118.2921, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)
  
})


test_that("fisher() works correctly with effective number of tests.", {
  
  res_nyh <- fisher(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj <- fisher(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- fisher(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- fisher(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  
  expect_equivalent(as.numeric(res_nyh$p), 5.26e-08, tolerance = p_tol)
  expect_equivalent(as.numeric(res_nyh$statistic), 113.1489, tolerance = stat_tol)

  expect_equivalent(as.numeric(res_lj$p), 1.00e-07, tolerance = p_tol)
  expect_equivalent(as.numeric(res_lj$statistic), 108.0058, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_gao$p), 2.75e-08, tolerance = p_tol)
  expect_equivalent(as.numeric(res_gao$statistic), 118.2921, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_gal$p), 1.92e-07, tolerance = p_tol)
  expect_equivalent(as.numeric(res_gal$statistic), 102.8627, tolerance = stat_tol)
  
})


test_that("fisher() works correctly with empirically-derived null distributions.", {
  
  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld)
  
  expect_equivalent(as.numeric(res$p), 0.00210, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 118.292, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0013, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0032, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)
  
  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)
  
  expect_equivalent(as.numeric(res$p), 0.00203, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 118.292, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0018, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0023, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)
  
  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)
  
  expect_equivalent(as.numeric(res$p), 0.00204, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 118.292, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0020, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0021, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))
  
  expect_equivalent(as.numeric(res$p), 0.00240, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 118.292, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0021, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0027, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)
  
})


test_that("fisher() works correctly under multivariate theory.", {
  
  res1 <- fisher(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 1))
  
  expect_equivalent(as.numeric(res1$p), 3.12e-05, tolerance = p_tol)
  expect_equivalent(as.numeric(res1$statistic), 62.90901, tolerance = stat_tol)
  expect_equivalent(attributes(res1$statistic)$df, 24.46, tolerance = df_tol)
  
  res2 <- fisher(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 2))
  
  expect_equivalent(as.numeric(res2$p), 0.000764, tolerance = p_tol)
  expect_equivalent(as.numeric(res2$statistic), 38.3482, tolerance = stat_tol)
  expect_equivalent(attributes(res2$statistic)$df, 14.91, tolerance = df_tol)
  
})
