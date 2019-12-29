
### library(library); library(testthat); Sys.setenv(NOT_CRAN="true")

p_tol <- 1e-06  # user-defined tolerance (for significant p-values)
stat_tol <- 1e-04  # user-defined tolerance (for test statistics and ci bounds)
df_tol <- 1e-02  # user-defined tolerance (for df)

context("Checking bonferroni() function")

test_that("bonferroni() works correctly under independence.", {
  
  res <- bonferroni(grid2ip.p)
  
  expect_equivalent(as.numeric(res$p), 0.068471, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 0.0030, tolerance = stat_tol)
  
})


test_that("bonferroni() works correctly with effective number of tests.", {
  
  res_nyh <- bonferroni(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj <- bonferroni(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- bonferroni(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- bonferroni(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  
  expect_equivalent(as.numeric(res_nyh$p), 0.065494, tolerance = p_tol)
  expect_equivalent(as.numeric(res_nyh$statistic), 0.0030, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_lj$p), 0.062517, tolerance = p_tol)
  expect_equivalent(as.numeric(res_lj$statistic), 0.0030, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_gao$p), 0.068471, tolerance = p_tol)
  expect_equivalent(as.numeric(res_gao$statistic), 0.0030, tolerance = stat_tol)
  
  expect_equivalent(as.numeric(res_gal$p), 0.05954, tolerance = p_tol)
  expect_equivalent(as.numeric(res_gal$statistic), 0.0030, tolerance = stat_tol)
  
})


test_that("bonferroni() works correctly with empirically-derived null distributions.", {
  
  set.seed(1234)
  res <- bonferroni(grid2ip.p, adjust = "empirical", R = grid2ip.ld)
  
  expect_equivalent(as.numeric(res$p), 0.049395, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 0.0030, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0452, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0538, tolerance = stat_tol)
  
  set.seed(1234)
  res <- bonferroni(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)
  
  expect_equivalent(as.numeric(res$p), 0.05245, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 0.0030, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0511, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0538, tolerance = stat_tol)
  
  set.seed(1234)
  res <- bonferroni(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)
  
  expect_equivalent(as.numeric(res$p), 0.051378, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 0.0030, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0509, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0518, tolerance = stat_tol)
  
  set.seed(1234)
  res <- bonferroni(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))
  
  expect_equivalent(as.numeric(res$p), 0.048395, tolerance = p_tol)
  expect_equivalent(as.numeric(res$statistic), 0.0030, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[1]), 0.0443, tolerance = stat_tol)
  expect_equivalent(as.numeric(res$ci[2]), 0.0528, tolerance = stat_tol)
  
})

