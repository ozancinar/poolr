
binotest <- function(p, adjust = "none", pca.method = NULL, R = NULL, alpha = 0.05, size = 10000, seed = NULL, ...) {
   if(adjust == "none") {
      k <- length(p)
      m <- sum(p <= alpha)
      probs <- dbinom(0:k, k, alpha)
      testStat <- dbinom(m, k, alpha)
      pooled.p <- sum(probs[probs <= testStat])
   } else if(adjust == "m.eff") {
      k <- length(p)
      m <- sum(p <= alpha)
      if(is.numeric(pca.method)) {
         eff <- pca.method 
      } else {
         eff <- meff(R = R, method = pca.method)
      }
      probs <- dbinom(0:eff, eff, alpha)
      testStat <- dbinom(round(m * eff / k), eff, alpha)
      pooled.p <- sum(probs[probs <= testStat])
   } else if (adjust == "empirical") {
      k <- length(p)
      m <- sum(p <= alpha)
      probs <- dbinom(0:k, k, alpha)
      testStat <- dbinom(m, k, alpha)
      tmp.p <- sum(probs[probs <= testStat])
      
      method <- "binotest"
      
      tmp <- list(...)
      if(is.null(tmp$emp.dis)) {
         emp.dist <- empirical(p = p, R = R, method = method, size = size, seed = seed)
      } else {
         emp.dist <- tmp$emp.dist
      }
      
      pooled.p <- sum(emp.dist <= testStat) / length(emp.dist)
   }
   
   res <- list(p = pooled.p, testStat = testStat, adjust = paste0(adjust, " ", pca.method))
   class(res) <- "combP"
   return(res)
}
