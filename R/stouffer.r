
stouffer <- function(p, adjust = "none", pca.method = NULL, R = NULL, size = 10000, seed = NULL, ...) {
   if(adjust == "none") {
      k <- length(p)
      testStat <- sum(qnorm(p)) / sqrt(k)
      pooled.p <- pnorm(testStat, lower.tail = FALSE)
   } else if(adjust == "m.eff") {
      k <- length(p)
      if(is.numeric(pca.method)) {
         eff <- pca.method 
      } else {
         eff <- meff(x = R, method = pca.method)
      }
      testStat <- sum(qnorm(p)) * sqrt(eff / k) / sqrt(k)
      pooled.p <- pnorm(testStat, lower.tail = FALSE)
   } else if (adjust == "empirical") {
      k <- length(p)
      testStat <- sum(qnorm(p)) / sqrt(k)
      tmp.p <- pnorm(testStat, lower.tail = FALSE)
      
      method <- "stouffer"
      
      tmp <- list(...)
      if(is.null(tmp$emp.dis)) {
         emp.dist <- empirical(p = p, R = R, method = method, size = size, seed = seed)
      } else {
         emp.dist <- tmp$emp.dist
      }
      
      pooled.p <- sum(emp.dist >= testStat) / length(emp.dist)
   } else if (adjust == "general") {
      k <- length(p)
      testStat <- sum(qnorm(p)) / sqrt(sum(R))
      pooled.p <- pnorm(testStat, lower.tail = FALSE)
   }
   
   res <- list(p = pooled.p, testStat = testStat, adjust = paste0(adjust, " ", pca.method))
   class(res) <- "combP"
   return(res)
}
