
empirical <- function(p, R = "ind", method, size = 10000, seed = "default", ...) {
   
   k <- length(p)
   
   if(R == "ind") {
      R <- diag(1, k)
   }
   
   if(seed != "default") {
      set.seed(seed)
   }
   
   if(method == "minp") {
      
      z <- mvrnorm(size, mu = rep(0, k), Sigma = R)
      pVals <- pnorm(z, lower.tail = FALSE)
      emp <- apply(pVals, 1, minp)
      
   } else if(method == "binotest") {
      
      z <- mvrnorm(size, mu = rep(0, k), Sigma = R)
      pVals <- pnorm(z, lower.tail = FALSE)
      emp <- apply(pVals, 1, binotest)
      
   } else if(method == "fisher") {
      
      z <- mvrnorm(size, mu = rep(0, k), Sigma = R)
      pVals <- pnorm(z, lower.tail = FALSE)
      emp <- apply(pVals, 1, fisher)
      
   } else if(method == "stouffer") {
      
      z <- mvrnorm(size, mu = rep(0, k), Sigma = R)
      pVals <- pnorm(z, lower.tail = FALSE)
      emp <- apply(pVals, 1, stouffer)
   }
   
   return(emp)
}