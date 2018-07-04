tippett <- function(p, adjust = "none", pca.method = NULL, R, size = 10000, seed = NULL, type = 2, ...) {

   k <- length(p)

   if (adjust == "none") {

      testStat <- 1 - (1 - min(p))^k
      pooled.p <- testStat
      adjust <- "none"

   }

   if (adjust == "m.eff") {

      if (is.numeric(pca.method)) {
         eff <- pca.method
         adjust <- paste0(pca.method, " (user defined)")
      } else {
         eff <- meff(x = R, method = pca.method)
         adjust <- paste0("meff (", pca.method, ")")
      }
      testStat <- 1 - (1 - min(p))^eff
      pooled.p <- testStat

   }

   if (adjust == "empirical") {

      tmp.p <- 1 - (1 - min(p))^k
      testStat <- 1 - (1 - min(p))^k
      method <- "tippett"

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(p = p, R = R, method = method, type = type, size = size, seed = seed)
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
