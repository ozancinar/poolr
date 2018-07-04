binotest <- function(p, adjust = "none", pca.method = NULL, R = NULL, alpha = 0.05, size = 10000, seed = NULL, type = 2, ...) {

   k <- length(p)
   r <- sum(p <= alpha)

   if (adjust == "none") {

      testStat <- dbinom(r, k, alpha)
      pooled.p <- sum(dbinom(r:k, k, alpha))
      adjust <- "none"

   } else if (adjust == "m.eff") {

      if (is.numeric(pca.method)) {
         eff <- pca.method
         adjust <- paste0(pca.method, " (user defined)")
      } else {
         eff <- meff(x = R, method = pca.method)
         adjust <- paste0("meff (", pca.method, ")")
      }
      testStat <- dbinom(round(r * eff / k), eff, alpha)
      pooled.p <- sum(dbinom(round(r * eff / k):eff, eff, alpha))

   } else if (adjust == "empirical") {

      probs <- dbinom(0:k, k, alpha)
      testStat <- dbinom(r, k, alpha)
      tmp.p <- sum(probs[probs <= testStat])

      method <- "binotest"

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(p = p, R = R, method = method, type = type, size = size, seed = seed)
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
