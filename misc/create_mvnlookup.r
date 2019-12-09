############################################################################

# Code to create the mvnlookup.rdata file.

# note: 15865 secs for "pracma" (with n=1000) on 'psysim' using 18 cores
#        2853 secs for "cubature" on 'psysim' using 18 cores

############################################################################

library(parallel)

cl <- makePSOCKcluster(18)

# vector of rho values for which we obtain the covariances
rhos <- seq(-0.99, 1, by = 0.01)

# choose method for the numerical integration
method <- "pracma"
#method <- "cubature"

# number of nodes for pracma::gaussLegendre()
n <- 500

# limits and tolerance for cubature::adaptIntegrate()
lims <- 5
tol  <- 1e-07

# value to round final results to
rnd <- 4

# load required packages in workers
invisible(clusterEvalQ(cl, {
   library(cubature)
   library(pracma)
   library(mvtnorm)
}))

# export 'n', 'lims', and 'tol'
clusterExport(cl, c("n", "lims", "tol"))

# set up data frame for storing results
mvnlookup <- data.frame(rhos = rhos)

time.start <- proc.time()

############################################################################

if (method == "pracma") {
   doint <- function(fun, xa = -5, xb = 5, ya = -5, yb = 5, n = 32, rho) {
      !any(is.numeric(xa), length(xa) == 1, is.numeric(ya), length(ya) == 1,
           is.numeric(xb), length(xb) == 1, is.numeric(yb), length(yb) == 1)
      cx <- gaussLegendre(n, xa, xb)
      x <- cx$x
      wx <- cx$w
      cy <- gaussLegendre(n, ya, yb)
      y <- cy$x
      wy <- cy$w
      mgrid <- meshgrid(x, y)
      Z <- matrix(NA, n, n)
      for(a in 1:n) {
         for(b in 1:n) {
            Z[a, b] <- fun(c(mgrid$X[a, b], mgrid$Y[a, b]), rho)
         }
      }
      Q <- wx %*% Z %*% as.matrix(wy)
      Q <- as.numeric(Q)
      return(Q)
   }
}

if (method == "cubature") {
   doint <- function(rho, tol, lims) {
      # adaptIntegrate() gets stuck in some cases when rho = 0, so skips this (we know that cov = 0 then anyway)
      if (rho == 0) {
         cov <- NA
      } else {
         # for 'z_2', using c(-lims,lims) leads to -Inf; make upper lims slightly larger to avoid this
         cov <- adaptIntegrate(intfun, tol=tol, lowerLimit=c(-lims,-lims), upperLimit=c(lims+.01,lims+.01), rho=rho)$integral
      }
      return(cov)
   }
}

clusterExport(cl, "doint")

############################################################################

# Cov(-2ln(p_i), -2ln(p_j)) for one-sided p-values

intfun <- function(xy, rho) {
   fx <- -2 * pnorm(xy[1], log.p = TRUE)
   fy <- -2 * pnorm(xy[2], log.p = TRUE)
   fx * fy * dmvnorm(xy, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
}

clusterExport(cl, "intfun")

if (method == "pracma")
   covs <- parSapplyLB(cl, rhos, function(r) doint(intfun, n=n, rho=r)) - 4
if (method == "cubature")
   covs <- parSapplyLB(cl, rhos, function(r) doint(rho=r, tol=tol, lims=lims)) - 4

covs[rhos == 1] <- 4
mvnlookup$m2lp_1 <- covs

############################################################################

# Cov(-2ln(p_i), -2ln(p_j)) for two-sided p-values

intfun <- function(xy,rho) {
   fx <- -2 * (log(2) + pnorm(abs(xy[1]), log.p = TRUE, lower.tail = FALSE))
   fy <- -2 * (log(2) + pnorm(abs(xy[2]), log.p = TRUE, lower.tail = FALSE))
   fx * fy * dmvnorm(xy, mean=c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
}

clusterExport(cl, "intfun")

if (method == "pracma")
   covs <- parSapplyLB(cl, rhos, function(r) doint(intfun, n=n, rho=r)) - 4
if (method == "cubature")
   covs <- parSapplyLB(cl, rhos, function(r) doint(rho=r, tol=tol, lims=lims)) - 4

covs[rhos == 1] <- 4
mvnlookup$m2lp_2 <- covs

############################################################################

# Cov(z_i, z_j) for one-sided p-values

intfun <- function(xy, rho) {
   fx <- qnorm(pnorm(xy[1]), lower.tail = FALSE)
   fy <- qnorm(pnorm(xy[2]), lower.tail = FALSE)
   fx * fy * dmvnorm(xy, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
}

clusterExport(cl, "intfun")

if (method == "pracma")
   covs <- parSapplyLB(cl, rhos, function(r) doint(intfun, n=n, rho=r))
if (method == "cubature")
   covs <- parSapplyLB(cl, rhos, function(r) doint(rho=r, tol=tol, lims=lims))

covs[rhos == 1] <- 1
mvnlookup$z_1 <- covs

############################################################################

# Cov(z_i, z_j) for two-sided p-values

intfun <- function(xy, rho) {
   fx <- qnorm(2 * pnorm(abs(xy[1]), lower.tail = FALSE), lower.tail = FALSE)
   fy <- qnorm(2 * pnorm(abs(xy[2]), lower.tail = FALSE), lower.tail = FALSE)
   fx * fy * dmvnorm(xy, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
}

clusterExport(cl, "intfun")

if (method == "pracma")
   covs <- parSapplyLB(cl, rhos, function(r) doint(intfun, n=n, rho=r))
if (method == "cubature")
   covs <- parSapplyLB(cl, rhos, function(r) doint(rho=r, tol=tol, lims=lims))

covs[rhos == 1] <- 1
mvnlookup$z_2 <- covs

############################################################################

# Cov(X^2_i, X^2_j) for one-sided p-values

intfun <- function(xy, rho) {
   fx <- qchisq(pnorm(xy[1]), df = 1, lower.tail = FALSE)
   fy <- qchisq(pnorm(xy[2]), df = 1, lower.tail = FALSE)
   fx * fy * dmvnorm(xy, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
}

clusterExport(cl, "intfun")

if (method == "pracma")
   covs <- parSapplyLB(cl, rhos, function(r) doint(intfun, n=n, rho=r)) - 1
if (method == "cubature")
   covs <- parSapplyLB(cl, rhos, function(r) doint(rho=r, tol=tol, lims=lims)) - 1

covs[rhos == 1] <- 2
mvnlookup$chisq1_1 <- covs

############################################################################

intfun <- function(xy, rho) {
   fx <- qchisq(2 * pnorm(abs(xy[1]), lower.tail = FALSE), df = 1, lower.tail = FALSE)
   fy <- qchisq(2 * pnorm(abs(xy[2]), lower.tail = FALSE), df = 1, lower.tail = FALSE)
   fx * fy * dmvnorm(xy, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
}

clusterExport(cl, "intfun")

if (method == "pracma")
   covs <- parSapplyLB(cl, rhos, function(r) doint(intfun, n=n, rho=r)) - 1
if (method == "cubature")
   covs <- parSapplyLB(cl, rhos, function(r) doint(rho=r, tol=tol, lims=lims)) - 1

covs[rhos == 1] <- 2
mvnlookup$chisq1_2 <- covs

############################################################################

# Cov(p_i, p_j) for one-sided p-values

intfun <- function(xy, rho) {
   fx <- pnorm(xy[1])
   fy <- pnorm(xy[2])
   fx * fy * dmvnorm(xy, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
}

clusterExport(cl, "intfun")

if (method == "pracma")
   covs <- parSapplyLB(cl, rhos, function(r) doint(intfun, n=n, rho=r)) - 1/4
if (method == "cubature")
   covs <- parSapplyLB(cl, rhos, function(r) doint(rho=r, tol=tol, lims=lims)) - 1/4

covs[rhos == 1] <- 1/12
mvnlookup$p_1 <- covs

############################################################################

# Cov(p_i, p_j) for two-sided p-values

intfun <- function(xy, rho) {
   fx <- 2 * pnorm(abs(xy[1]), lower.tail = FALSE)
   fy <- 2 * pnorm(abs(xy[2]), lower.tail = FALSE)
   fx * fy * dmvnorm(xy, mean = c(0, 0), sigma = matrix(c(1, rho, rho, 1), nrow = 2))
}

clusterExport(cl, "intfun")

if (method == "pracma")
   covs <- parSapplyLB(cl, rhos, function(r) doint(intfun, n=n, rho=r)) - 1/4
if (method == "cubature")
   covs <- parSapplyLB(cl, rhos, function(r) doint(rho=r, tol=tol, lims=lims)) - 1/4

covs[rhos == 1] <- 1/12
mvnlookup$p_2 <- covs

############################################################################

stopCluster(cl)

time.end <- proc.time()
secs <- unname(time.end - time.start)[3]

# total running time in seconds and minutes
print(round(secs))
print(round(secs / 60, 1))

# set all covs to 0 for rho = 0
mvnlookup[rhos == 0,] <- 0

# round results
mvnlookup <- round(mvnlookup, rnd)

# save results
save(mvnlookup, file = "../data/mvnlookup.rdata")
#save(mvnlookup, file = paste0("mvnlookup_", method, ".rdata"))

############################################################################

if (F) {

   par(mfrow=c(4, 2))

   for (i in 1:8) {
      print(names(mvnlookup)[i+1])
      x <- mvnlookup$rhos
      y <- mvnlookup[[i + 1]]
      X <- poly(x, raw=TRUE, simple=TRUE, degree=8)
      colnames(X) <- paste0(".", 1:ncol(X))
      weights <- c(rep(1, length(x) - 1), 1000)
      res <- lm(y ~ 0 + X, weights=weights)
      print(summary(res))
      print(round(coef(res), 4))

      if (F) {
         plot(x, y, pch=19, cex=0.2, main=names(mvnlookup)[i+1], xlab="rho", ylab="cov")
         abline(h=0, lty="dotted")
         abline(v=0, lty="dotted")
         lines(x, fitted(res))
      } else {
         abslim <- .01
         plot(x, resid(res), pch=19, cex=0.2, main=names(mvnlookup)[i+1], type="l", ylim=c(-abslim,abslim), xlab="rho", ylab="resid")
         abline(h=0, lty="dotted")
         abline(v=0, lty="dotted")
      }

   }

}

############################################################################
