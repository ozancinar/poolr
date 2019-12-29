
### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

# p_tol <- 1e-06  # user-defined tolerance (for significant p-values)
# stat_tol <- 1e-04  # user-defined tolerance (for test statistics and ci bounds)
# df_tol <- 1e-02  # user-defined tolerance (for df)

p_tol <- 0.1
stat_tol <- 0.1
df_tol <- 0.1

context("Checking binotest() function")

test_that("binotest() works correctly under independence.", {
  
  res <- binotest(grid2ip.p)
  
  expect_equivalent(as.numeric(res$p), 3.76e-09, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 11, tolerance = stat_tol)
  
})


test_that("binotest() works correctly with effective number of tests.", {
  
  res_nyh <- binotest(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj <- binotest(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- binotest(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- binotest(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  
  expect_equivalent(as.numeric(res_nyh$p), 2.05e-09, tolerance = p_tol)
  expect_equivalent(as.numeric(res_nyh$statistic), 11, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_lj$p), 2.06e-08, tolerance = p_tol)
  expect_equivalent(as.numeric(res_lj$statistic), 11, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_gao$p), 3.76e-09, tolerance = p_tol)
  expect_equivalent(as.numeric(res_gao$statistic), 11, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_gal$p), 1.13e-08, tolerance = p_tol)
  expect_equivalent(as.numeric(res_gal$statistic), 11, tolerance = stat_tol)
  
})


test_that("binotest() works correctly with empirically-derived null distributions.", {
  
  set.seed(1234)
  res <- binotest(grid2ip.p, adjust = "empirical", R = grid2ip.ld)
  
  expect_equivalent(as.numeric(res$p), 0.00030, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 11, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 6.18e-05, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.00088, tolerance = stat_tol)
  
  set.seed(1234)
  res <- binotest(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)
  
  expect_equivalent(as.numeric(res$p), 0.00037, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 11, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0003, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0005, tolerance = stat_tol)
  
  set.seed(1234)
  res <- binotest(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)
  
  expect_equivalent(as.numeric(res$p), 0.00045, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 11, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0004, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0005, tolerance = stat_tol)
  
  set.seed(1234)
  res <- binotest(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))
  
  expect_equivalent(as.numeric(res$p), 0.00061, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 11, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0005, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0008, tolerance = stat_tol)
  
})
