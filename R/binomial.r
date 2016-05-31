
binomial <- function(p, adjust = "none", pca.method = NULL, R = NULL, ...) {
   if(adjust == "none") {
      k <- length(p)
      m <- sum(p <= 0.05)
      probs <- dbinom(0:k, k, 0.05)
      prob <- dbinom(m, k, 0.05)
      pooled.p <- sum(probs[probs <= prob])
   } else if(adjust == "m.eff") {
      k <- length(p)
      m <- sum(p <= 0.05)
      eff <- m.eff(R = R, method = pca.method)
      probs <- dbinom(0:eff, eff, 0.05)
      prob <- dbinom(m, eff, 0.05)
      pooled.p <- sum(probs[probs <= prob])
   }
   res <- list(p = pooled.p, adjust = paste0(adjust, " ", pca.method))
   return(res)
}
