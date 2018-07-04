invchisq <- function(p, adjust = "none", pca.method = NULL, R, size = 10000, seed = NULL, type = 2, ...) {

   k <- length(p)

   if (adjust == "none") {

      testStat <- sum(qchisq(p, df = 1, lower.tail = FALSE))
      pooled.p <- pchisq(testStat, df = k, lower.tail=FALSE)
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
      testStat <- sum(qchisq(p, df = 1, lower.tail = FALSE)) * (eff / k)
      pooled.p <- pchisq(testStat, df = eff, lower.tail = FALSE)

  }

  if (adjust == "empirical") {

      testStat <- sum(qchisq(p, df = 1, lower.tail = FALSE))
      tmp.p <- pchisq(testStat, df = k, lower.tail=FALSE)
      method <- "invchisq"

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(p = p, R = R, method = method, type = type, size = size, seed = seed)
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
