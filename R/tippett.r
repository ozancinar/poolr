tippett <- function(p, adjust = "none", pca.method, R, size = 10000, seed, type = 2, ...) {

   k <- length(p)

   testStat <- 1 - (1 - min(p))^k

   if (adjust == "none") {

      pooled.p <- testStat
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

      testStat <- 1 - (1 - min(p))^m
      pooled.p <- testStat

   }

   if (adjust == "empirical") {

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = "tippett", type = type, size = size, seed = seed)
      } else {
         emp.dist <- tmp$emp.dist
      }

      pooled.p <- sum(emp.dist <= testStat) / length(emp.dist)
      adjust <- "empirical"

   }

   res <- list(p = pooled.p, testStat = testStat, adjust = adjust)
   class(res) <- "combP"
   return(res)

}
