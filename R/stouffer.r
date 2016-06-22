
stouffer <- function(p, adjust = "none", pca.method = NULL, R = NULL, size = 10000, seed = NULL, ...) {
   if(adjust == "none") {
      k <- length(p)
      pooled.p <- 2 * pnorm(abs(sum(qnorm(p)) / sqrt(k)), lower.tail = FALSE)
   } else if(adjust == "m.eff") {
      k <- length(p)
      eff <- meff(R = R, method = pca.method)
      pooled.p <- 2 * pnorm(abs(sum(qnorm(p)) / sqrt(eff)), lower.tail = FALSE)
   } else if (adjust == "empirical") {
      k <- length(p)
      tmp.p <- 2 * pnorm(abs(sum(qnorm(p)) / sqrt(k)), lower.tail = FALSE)
      
      method <- "stouffer"
      
      emp.dist <- empirical(p = p, R = R, method = method, size = size, seed = seed)
      pooled.p <- sum(emp.dist <= tmp.p) / length(emp.dist)
   }
   
   res <- list(p = pooled.p, adjust = paste0(adjust, " ", pca.method))
   return(res)
}
