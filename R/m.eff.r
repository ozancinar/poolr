
m.eff <- function(R, method = NULL, ...) {
  
  # number of p-values
  k <- nrow(R)
  
  # dimension checks
  if(!isSymmetric(R)) {
    stop("R is not symmetric.")
  }
  
  if(nrow(R) != k) {
    stop("Dimensions of R do not match with the length of p.")
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
     m.eff <- 1 + (k - 1) * (1 - var(evs) / k)
     m.eff <- floor(m.eff)
  } else if(method == "li.ji") {
     ### effective number of tests (based on Li & Ji, 2005)
     m.eff <- sum(ifelse(abs.evs >= 1, 1, 0) + (abs.evs - floor(abs.evs)))
     m.eff <- floor(m.eff)
  } else if(method == "gao") {
     ### effective number of tests (based on Gao, 2008)
     m.eff <- which(cumsum(sort(abs.evs, decreasing = TRUE)) / sum(abs.evs) > 0.995)[1]
  } else if(method == "galwey") {
     ### effective number of tests (based on Galwey, 2009)
     evs[evs < 0] <- 0
     m.eff <- sum(sqrt(evs))^2 / sum(evs)
     m.eff <- floor(m.eff)
  }
  
  return(m.eff)
}
