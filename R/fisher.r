fisher <- function(p, adjust = "none", pca.method = NULL, R = NULL, ...) {
   if(adjust == "none") {
      k <- length(p)
      pooled.p <- pchisq(-2*sum(log(p)), df=2*k, lower.tail=FALSE)
   } else if(adjust == "m.eff") {
      k <- length(p)
      m.eff <- m.eff(p = p, R = R, method = pca.method)
      pooled.p <- pchisq(-2 * sum(log(p))/ m.eff, df = 2 * m.eff, lower.tail = FALSE)
   }
   res <- list(p = pooled.p, adjust = paste0(adjust, " ", pca.method))
   return(res)
}
