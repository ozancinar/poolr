
empirical <- function(p, R = NULL, method, size = 10000, seed = NULL, ...) {
   
   k <- length(p)
   
   if(is.null(R)) {
      R <- diag(1, k)
   }
   
   R <- nearPD(R)$mat
   
   if(!is.null(seed)) {
      set.seed(seed)
   }
   
   if(method == "minp") {
      
      z <- mvrnorm(size, mu = rep(0, k), Sigma = R)
      pVals <- pnorm(z, lower.tail = FALSE)
      emp <- apply(pVals, 1, function(x) {minp(x)$p})
      
   } else if(method == "binotest") {
      
      z <- mvrnorm(size, mu = rep(0, k), Sigma = R)
      pVals <- pnorm(z, lower.tail = FALSE)
      emp <- apply(pVals, 1, function(x) {binotest(x)$p})
      
   } else if(method == "fisher") {
      
      z <- mvrnorm(size, mu = rep(0, k), Sigma = R)
      pVals <- pnorm(z, lower.tail = FALSE)
      emp <- apply(pVals, 1, function(x) {fisher(x)$p})
      
   } else if(method == "stouffer") {
      
      z <- mvrnorm(size, mu = rep(0, k), Sigma = R)
      pVals <- pnorm(z, lower.tail = FALSE)
      emp <- apply(pVals, 1, function(x) {stouffer(x)$p})
   }
   
   return(emp)
}