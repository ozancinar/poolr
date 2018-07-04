stouffer <- function(p, adjust = "none", pca.method = NULL, R, size = 10000, seed, type = 2, ...) {

   k <- length(p)

   if (adjust == "none") {

      testStat <- sum(qnorm(p, lower.tail = FALSE)) / sqrt(k)
      pooled.p <- pnorm(testStat, lower.tail = FALSE)
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
      testStat <- sum(qnorm(p, lower.tail = FALSE)) * sqrt(eff / k) / sqrt(k)
      pooled.p <- pnorm(testStat, lower.tail = FALSE)

   }

   if (adjust == "general") {

      testStat <- sum(qnorm(p, lower.tail = FALSE)) / sqrt(sum(R))
      pooled.p <- pnorm(testStat, lower.tail = FALSE)
      adjust <- "generalized stouffer"

   }

   if (adjust == "empirical") {

      testStat <- sum(qnorm(p, lower.tail = FALSE)) / sqrt(k)
      method <- "stouffer"

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
