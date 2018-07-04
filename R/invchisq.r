invchisq <- function(p, adjust = "none", pca.method, R, size = 10000, seed, type = 2, ...) {

   k <- length(p)

   testStat <- sum(qchisq(p, df = 1, lower.tail = FALSE))

   if (adjust == "none") {

      pooled.p <- pchisq(testStat, df = k, lower.tail = FALSE)
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
      pooled.p <- pchisq(testStat, df = m, lower.tail = FALSE)

   }

   if (adjust == "empirical") {

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = "invchisq", type = type, size = size, seed = seed)
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
