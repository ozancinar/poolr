
minp <- function(p, adjust = "none", pca.method = NULL, R = NULL, ...) {
   if (adjust == "none") {
      k <- length(p)
      pooled.p <- min(p) * k
   } else if (adjust == "m.eff") {
      eff <- meff(R = R, method = pca.method)
      pooled.p <- min(p) * eff
   }
   res <- list(p = pooled.p, adjust = paste0(adjust, " "))
   return(res)
}