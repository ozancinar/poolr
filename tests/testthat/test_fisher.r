### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")

source("tolerances.r")

context("Checking fisher() function")

test_that("fisher() works correctly under independence.", {

  res <- fisher(grid2ip.p)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 1.389547e-09, tolerance = p_tol)
  expect_equivalent(c(res$statistic), 127.4818, tolerance = stat_tol)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

})

test_that("fisher() works correctly with effective number of tests.", {

  res_nyh <- fisher(grid2ip.p, adjust = "nyholt", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_lj  <- fisher(grid2ip.p, adjust = "liji", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gao <- fisher(grid2ip.p, adjust = "gao", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_gal <- fisher(grid2ip.p, adjust = "galwey", R = mvnconv(grid2ip.ld, target = "p", cov2cor = TRUE))
  res_user <- fisher(grid2ip.p, m = 18)

  out <- capture.output(print(res_nyh))
  out <- capture.output(print(res_lj))
  out <- capture.output(print(res_gao))
  out <- capture.output(print(res_gal))
  out <- capture.output(print(res_user))

  expect_equivalent(c(res_nyh$p), 3.008722e-09, tolerance = p_tol)
  expect_equivalent(c(res_nyh$statistic), 121.9391, tolerance = stat_tol)

  expect_equivalent(c(res_lj$p), 6.52039e-09, tolerance = p_tol)
  expect_equivalent(c(res_lj$statistic), 116.3964, tolerance = stat_tol)

  expect_equivalent(c(res_gao$p), 1.389547e-09, tolerance = p_tol)
  expect_equivalent(c(res_gao$statistic), 127.4818, tolerance = stat_tol)

  expect_equivalent(c(res_gal$p), 1.414432e-08, tolerance = p_tol)
  expect_equivalent(c(res_gal$statistic), 110.8537, tolerance = stat_tol)

  expect_equivalent(c(res_user$p), 6.677494e-08, tolerance = p_tol)
  expect_equivalent(c(res_user$statistic), 99.76835, tolerance = stat_tol)

})

test_that("fisher() works correctly with empirically-derived null distributions.", {

  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.00049995, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 127.4818, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0001623517, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.001166328, tolerance = p_tol * emp_sca)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 100000)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.00101999, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 127.4818, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.000831754, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.001238063, tolerance = p_tol * emp_sca)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = 1000000, batchsize = 1000)
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.000972999, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 127.4818, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0009128416, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.001036076, tolerance = p_tol * emp_sca)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

  set.seed(1234)
  res <- fisher(grid2ip.p, adjust = "empirical", R = grid2ip.ld, size = c(1000, 10000, 100000), threshold = c(0.10, 0.01))
  out <- capture.output(print(res))

  expect_equivalent(c(res$p), 0.00104999, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$statistic), 127.4818, tolerance = stat_tol * emp_sca)
  expect_equivalent(c(res$ci[1]), 0.0008588652, tolerance = p_tol * emp_sca)
  expect_equivalent(c(res$ci[2]), 0.001270936, tolerance = p_tol * emp_sca)
  expect_equivalent(attributes(res$statistic)$df, 46, tolerance = df_tol)

})

test_that("fisher() works correctly under multivariate theory.", {

  res1 <- fisher(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 1))
  out <- capture.output(print(res1))

  expect_equivalent(c(res1$p), 6.063446e-06, tolerance = p_tol)
  expect_equivalent(c(res1$statistic), 67.80295, tolerance = stat_tol)
  expect_equivalent(attributes(res1$statistic)$df, 24.46574, tolerance = df_tol)

  res2 <- fisher(grid2ip.p, adjust = "generalized", R = mvnconv(grid2ip.ld, side = 2))
  out <- capture.output(print(res2))

  expect_equivalent(c(res2$p), 0.0002622587, tolerance = p_tol)
  expect_equivalent(c(res2$statistic), 41.55411, tolerance = stat_tol)
  expect_equivalent(attributes(res2$statistic)$df, 14.99421, tolerance = df_tol)

})
