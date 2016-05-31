
stouffer <- function(p, adjust = "none", pca.method = NULL, R = NULL, ...) {
   if(adjust == "none") {
      k <- length(p)
      pooled.p <- 2 * pnorm(abs(sum(qnorm(p)) / sqrt(k)), lower.tail = FALSE)
   } else if(adjust == "m.eff") {
      k <- length(p)
      eff <- m.eff(R = R, method = pca.method)
      pooled.p <- 2 * pnorm(abs(sum(qnorm(p)) / sqrt(eff)), lower.tail = FALSE)
   }
   res <- list(p = pooled.p, adjust = paste0(adjust, " ", pca.method))
   return(res)
}
