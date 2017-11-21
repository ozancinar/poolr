
binotest <- function(p, adjust = "none", pca.method = NULL, R = NULL, alpha = 0.05, size = 10000, seed = NULL, type = 2, ...) {
   if(adjust == "none") {
      k <- length(p)
      m <- sum(p <= alpha)
      testStat <- dbinom(m, k, alpha)
      pooled.p <- sum(dbinom(m:k, k, alpha))
      adjust <- "none"
   } else if(adjust == "m.eff") {
      k <- length(p)
      m <- sum(p <= alpha)
      if(is.numeric(pca.method)) {
         eff <- pca.method 
         adjust <- paste0(pca.method, " (user defined)")
      } else {
         eff <- meff(x = R, method = pca.method)
         adjust <- paste0("meff (", pca.method, ")")
      }
      testStat <- dbinom(round(m * eff / k), eff, alpha)
      pooled.p <- sum(dbinom(round(m * eff / k):eff, eff, alpha))
   } else if (adjust == "empirical") {
      k <- length(p)
      m <- sum(p <= alpha)
      probs <- dbinom(0:k, k, alpha)
      testStat <- dbinom(m, k, alpha)
      tmp.p <- sum(probs[probs <= testStat])
      
      method <- "binotest"
      
      tmp <- list(...)
      if(is.null(tmp$emp.dis)) {
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
