binotest <- function(p, adjust = "none", pca.method = NULL, R, alpha = 0.05, size = 10000, seed = NULL, type = 2, ...) {

   k <- length(p)
   r <- sum(p <= alpha)

   if (adjust == "none") {

      testStat <- dbinom(r, k, alpha)
      pooled.p <- sum(dbinom(r:k, k, alpha))
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
      testStat <- dbinom(round(r * eff / k), eff, alpha)
      pooled.p <- sum(dbinom(round(r * eff / k):eff, eff, alpha))

   }

   if (adjust == "empirical") {

      testStat <- dbinom(r, k, alpha)
      method <- "binotest"

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = method, type = type, size = size, seed = seed)
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
