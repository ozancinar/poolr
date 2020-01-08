### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")
library(poolr)

source("tolerances.r")

context("Checking meff() function")

test_that("meff() works correctly.", {

  m_nyh1 <- meff(R = mvnconv(grid2ip.ld, side = 1, target = "p", cov2cor = TRUE), method = "nyholt")
  m_lj1  <- meff(R = mvnconv(grid2ip.ld, side = 1, target = "p", cov2cor = TRUE), method = "liji")
  m_gao1 <- meff(R = mvnconv(grid2ip.ld, side = 1, target = "p", cov2cor = TRUE), method = "gao")
  m_gal1 <- meff(R = mvnconv(grid2ip.ld, side = 1, target = "p", cov2cor = TRUE), method = "galwey")

  expect_equivalent(m_nyh1, 21, tolerance = df_m)
  expect_equivalent(m_lj1,  15, tolerance = df_m)
  expect_equivalent(m_gao1, 20, tolerance = df_m)
  expect_equivalent(m_gal1, 14, tolerance = df_m)

  m_nyh2 <- meff(R = mvnconv(grid2ip.ld, side = 2, target = "p", cov2cor = TRUE), method = "nyholt")
  m_lj2  <- meff(R = mvnconv(grid2ip.ld, side = 2, target = "p", cov2cor = TRUE), method = "liji")
  m_gao2 <- meff(R = mvnconv(grid2ip.ld, side = 2, target = "p", cov2cor = TRUE), method = "gao")
  m_gal2 <- meff(R = mvnconv(grid2ip.ld, side = 2, target = "p", cov2cor = TRUE), method = "galwey")

  expect_equivalent(m_nyh2, 22, tolerance = df_m)
  expect_equivalent(m_lj2,  21, tolerance = df_m)
  expect_equivalent(m_gao2, 23, tolerance = df_m)
  expect_equivalent(m_gal2, 20, tolerance = df_m)

})
