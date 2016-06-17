
binotest <- function(p, adjust = "none", pca.method = NULL, R = NULL, alpha = 0.05, ...) {
   if(adjust == "none") {
      k <- length(p)
      m <- sum(p <= alpha)
      probs <- dbinom(0:k, k, alpha)
      prob <- dbinom(m, k, alpha)
      pooled.p <- sum(probs[probs <= prob])
   } else if(adjust == "m.eff") {
      k <- length(p)
      m <- sum(p <= alpha)
      eff <- meff(R = R, method = pca.method)
      probs <- dbinom(0:eff, eff, alpha)
      prob <- dbinom(round(m * eff / k), eff, alpha)
      pooled.p <- sum(probs[probs <= prob])
   }
   res <- list(p = pooled.p, adjust = paste0(adjust, " ", pca.method))
   return(res)
}
