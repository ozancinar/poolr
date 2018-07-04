fisher <- function(p, adjust = "none", pca.method = NULL, R, size = 10000, seed = NULL, type = 2, ...) {

   k <- length(p)

   if (adjust == "none") {

      testStat <- -2*sum(log(p))
      pooled.p <- pchisq(testStat, df=2*k, lower.tail = FALSE)
      adjust <- "none"

   }

   if (adjust == "m.eff") {

      if (is.numeric(pca.method)) {
         eff <- pca.method
         adjust <- paste0(pca.method, " (user defined)")
      } else {
         eff <- meff(R = R, method = pca.method)
         adjust <- paste0("meff (", pca.method, ")")
      }
      testStat <- -2 * sum(log(p)) * (eff / k)
      pooled.p <- pchisq(-2 * sum(log(p)) * (eff / k), df = 2 * eff, lower.tail = FALSE)

   }

   if (adjust == "brown") {

      tmp <- list(...)
      if (is.null(tmp$brownCov)) {
         covs <- brown(R)
      } else {
         covs <- tmp$brownCov
      }

      chi2val <- -2 * sum(log(p))
      expx2 <- 2 * k
      #varx2 <- 4 * k + 2 * sum(covs[upper.tri(covs)])
      varx2 <- sum(covs)
      fval <- 2 * expx2^2 / varx2
      cval <- varx2 / (2 * expx2)

      testStat <- chi2val/cval
      pooled.p <- pchisq(chi2val/cval, df=fval, lower.tail = FALSE)
      adjust <- "brown"

   }

   if (adjust == "empirical") {

      testStat <- -2*sum(log(p))
      method <- "fisher"

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = method, type = type, size = size, seed = seed)
      } else {
         emp.dist <- tmp$emp.dist
      }

      pooled.p <- sum(emp.dist >= testStat) / length(emp.dist)
      adjust <- "empirical"

   }

   res <- list(p = pooled.p, testStat = testStat, adjust = adjust)
   class(res) <- "combP"
   return(res)

}
