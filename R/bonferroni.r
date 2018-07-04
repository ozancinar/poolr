bonferroni <- function(p, adjust = "none", pca.method = NULL, R, size = 10000, seed = NULL, type = 2, ...) {

   k <- length(p)

   if (adjust == "none") {

      testStat <- min(p) * k
      pooled.p <- testStat
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
      testStat <- min(p) * eff
      pooled.p <- min(p) * eff

   }

   if (adjust == "empirical") {

      testStat <- min(p) * k
      method <- "bonferroni"

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = method, type = type, size = size, seed = seed)
      } else {
         emp.dist <- tmp$emp.dist
      }

      pooled.p <- sum(emp.dist <= testStat) / length(emp.dist)
      adjust <- "empirical"

   }

   if (pooled.p > 1)
      pooled.p <- 1

   res <- list(p = pooled.p, testStat = testStat, adjust = adjust)
   class(res) <- "combP"
   return(res)

}
