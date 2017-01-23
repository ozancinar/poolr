
meff <- function(R, method = NULL, ...) {
  
  # number of p-values
  k <- nrow(R)
  
  # dimension checks
  if(!isSymmetric(R)) {
    stop("R is not symmetric.")
  }
  
  if(!method %in% c("nyholt", "li.ji", "gao", "galwey")) {
     stop("The method for PCA is not correct.")
  }
  
  # Converting the correlation matrix to positive-definite.
  R <- as.matrix(nearPD(R)$mat)
  
  ### get eigenvalues and absolute eigenvalues of R matrix
  evs <- eigen(R)$values
  abs.evs <- abs(evs)
  
  if(method == "nyholt") {
     ### effective number of tests (based on Nyholt, 2004)
     eff <- 1 + (k - 1) * (1 - var(evs) / k)
     eff <- floor(eff)
  } else if(method == "li.ji") {
     ### effective number of tests (based on Li & Ji, 2005)
     eff <- sum(ifelse(abs.evs >= 1, 1, 0) + (abs.evs - floor(abs.evs)))
     eff <- floor(eff)
  } else if(method == "gao") {
     ### effective number of tests (based on Gao, 2008)
     eff <- which(cumsum(sort(abs.evs, decreasing = TRUE)) / sum(abs.evs) > 0.995)[1]
  } else if(method == "galwey") {
     ### effective number of tests (based on Galwey, 2009)
     evs[evs < 0] <- 0
     eff <- sum(sqrt(evs))^2 / sum(evs)
     eff <- floor(eff)
  }
  
  return(eff)
}
