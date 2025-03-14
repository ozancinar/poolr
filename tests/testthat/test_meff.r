### library(poolr); library(testthat); Sys.setenv(NOT_CRAN="true")

source("tolerances.r")

context("Checking meff() function")

test_that("meff() works correctly.", {

  m_nyh1 <- meff(R = mvnconv(grid2ip.ld, side = 1, target = "p", cov2cor = TRUE), method = "nyholt")
  m_lj1  <- meff(R = mvnconv(grid2ip.ld, side = 1, target = "p", cov2cor = TRUE), method = "liji")
  m_gao1 <- meff(R = mvnconv(grid2ip.ld, side = 1, target = "p", cov2cor = TRUE), method = "gao")
  m_gal1 <- meff(R = mvnconv(grid2ip.ld, side = 1, target = "p", cov2cor = TRUE), method = "galwey")

  expect_equivalent(m_nyh1, 21, tolerance = m_tol)
  expect_equivalent(m_lj1,  15, tolerance = m_tol)
  expect_equivalent(m_gao1, 20, tolerance = m_tol)
  expect_equivalent(m_gal1, 14, tolerance = m_tol)

  m_nyh2 <- meff(R = mvnconv(grid2ip.ld, side = 2, target = "p", cov2cor = TRUE), method = "nyholt")
  m_lj2  <- meff(R = mvnconv(grid2ip.ld, side = 2, target = "p", cov2cor = TRUE), method = "liji")
  m_gao2 <- meff(R = mvnconv(grid2ip.ld, side = 2, target = "p", cov2cor = TRUE), method = "gao")
  m_gal2 <- meff(R = mvnconv(grid2ip.ld, side = 2, target = "p", cov2cor = TRUE), method = "galwey")

  expect_equivalent(m_nyh2, 22, tolerance = m_tol)
  expect_equivalent(m_lj2,  21, tolerance = m_tol)
  expect_equivalent(m_gao2, 23, tolerance = m_tol)
  expect_equivalent(m_gal2, 20, tolerance = m_tol)
  
  m_lj_eigen <- meff(eigen = eigen(mvnconv(grid2ip.ld, side = 2, target = "p", cov2cor = TRUE))$values, method = "liji")
  expect_equivalent(m_lj_eigen, 21, tolerance = m_tol)

})

test_that("The arguments of meff() are checked correctly.", {
  
  expect_error(meff(method = "liji"), "Argument 'R' must be specified.")
  expect_error(meff(eigen = c("a"), method = "liji"), "Argument 'eigen' must be a numeric vector.")
  expect_warning(meff(mvnconv(grid2ip.ld, side = 2, target = "p", cov2cor = TRUE), method = "gao", C = "a"))

})
