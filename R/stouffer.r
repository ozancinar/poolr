
stouffer <- function(p, adjust = "none", pca.method = NULL, R = NULL, size = 10000, seed = NULL, type = 2, ...) {
   if(adjust == "none") {
      k <- length(p)
      testStat <- sum(qnorm(p, lower.tail = FALSE)) / sqrt(k)
      pooled.p <- pnorm(testStat, lower.tail = FALSE)
      adjust <- "none"
   } else if(adjust == "m.eff") {
      k <- length(p)
      if(is.numeric(pca.method)) {
         eff <- pca.method 
         adjust <- paste0(pca.method, " (user defined)")
      } else {
         eff <- meff(x = R, method = pca.method)
         adjust <- paste0("meff (", pca.method, ")")
      }
      testStat <- sum(qnorm(p, lower.tail = FALSE)) * sqrt(eff / k) / sqrt(k)
      pooled.p <- pnorm(testStat, lower.tail = FALSE)
   } else if (adjust == "general") {
     k <- length(p)
     testStat <- sum(qnorm(p, lower.tail = FALSE)) / sqrt(sum(R))
     pooled.p <- pnorm(testStat, lower.tail = FALSE)
     adjust <- "generalized stouffer"
   } else if (adjust == "empirical") {
      k <- length(p)
      testStat <- sum(qnorm(p, lower.tail = FALSE)) / sqrt(k)
      tmp.p <- pnorm(testStat, lower.tail = FALSE)
      
      method <- "stouffer"
      
      tmp <- list(...)
      if(is.null(tmp$emp.dis)) {
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
