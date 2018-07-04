fisher <- function(p, adjust = "none", pca.method, R, size = 10000, seed, type = 2, ...) {

   k <- length(p)

   testStat <- -2 * sum(log(p))

   if (adjust == "none") {

      pooled.p <- pchisq(testStat, df = 2*k, lower.tail = FALSE)
      adjust <- "none"

   }

   if (adjust == "m.eff") {

      if (is.numeric(pca.method)) {
         m <- pca.method
         adjust <- paste0(pca.method, " (user defined)")
      } else {
         m <- meff(R = R, method = pca.method)
         adjust <- paste0("meff (", pca.method, ")")
      }

      testStat <- testStat * (m / k)
      pooled.p <- pchisq(testStat * (m / k), df = 2 * m, lower.tail = FALSE)

   }

   if (adjust == "brown") {

      tmp <- list(...)
      if (is.null(tmp$brownCov)) {
         covs <- brown(R)
      } else {
         covs <- tmp$brownCov
      }

      expx2 <- 2 * k
      varx2 <- sum(covs)
      fval <- 2 * expx2^2 / varx2
      cval <- varx2 / (2 * expx2)

      testStat <- testStat/cval
      pooled.p <- pchisq(testStat, df = fval, lower.tail = FALSE)
      adjust <- "brown"

   }

   if (adjust == "empirical") {

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = "fisher", type = type, size = size, seed = seed)
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
