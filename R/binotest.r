binotest <- function(p, adjust = "none", pca.method, R, alpha = 0.05, size = 10000, seed, type = 2, ...) {

   k <- length(p)
   r <- sum(p <= alpha)

   testStat <- dbinom(r, k, alpha)

   if (adjust == "none") {

      pooled.p <- sum(dbinom(r:k, k, alpha))
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

      testStat <- dbinom(round(r * m / k), m, alpha)
      pooled.p <- sum(dbinom(round(r * m / k):m, m, alpha))

   }

   if (adjust == "empirical") {

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = "binotest", type = type, size = size, seed = seed)
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
