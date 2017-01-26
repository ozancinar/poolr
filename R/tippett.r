
tippett <- function(p, adjust = "none", pca.method = NULL, R = NULL, size = 10000, seed = NULL, ...) {
   if (adjust == "none") {
      k <- length(p)
      testStat <- 1 - (1 - min(p))^k
      pooled.p <- 1 - (1 - min(p))^k
   } else if (adjust == "m.eff") {
      if(is.numeric(pca.method)) {
         eff <- pca.method 
      } else {
         eff <- meff(R = R, method = pca.method)
      }
      testStat <- 1 - (1 - min(p))^eff
      pooled.p <- 1 - (1 - min(p))^eff
   } else if (adjust == "empirical") {
      k <- length(p)
      tmp.p <- 1 - (1 - min(p))^k
      testStat <- 1 - (1 - min(p))^k
      method <- "tippett"
      
      tmp <- list(...)
      if(is.null(tmp$emp.dis)) {
         emp.dist <- empirical(p = p, R = R, method = method, size = size, seed = seed)
      } else {
         emp.dist <- tmp$emp.dist
      }
      
      pooled.p <- sum(emp.dist <= testStat) / length(emp.dist)
   }
   
   if(pooled.p > 1) {pooled.p <- 1}
   res <- list(p = pooled.p, testStat = testStat, adjust = paste0(adjust, " ", pca.method))
   class(res) <- "combP"
   return(res)
}