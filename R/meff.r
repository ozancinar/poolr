meff <- function(R, eigen = FALSE, method, ...) {
  
  # match method argument
  method <- match.arg(method, c("nyholt", "liji", "gao", "galwey"))
  
  if (eigen) {
    
    # can pass eigenvalues directly to function if eigen is TRUE
    
    if (!class(R) %in% c("numeric", "integer"))
      stop("eigenvalues are not numeric or integer.")
    
    evs <- R
    abs.evs <- abs(evs)
    
  } else {
    
    # check that R is symmetric
    if (!isSymmetric(R))
      stop("R is not symmetric.")
    
    # ensure that the correlation matrix is positive semi-definite
    #R <- as.matrix(nearPD(R)$mat)
    
    # get eigenvalues and absolute eigenvalues of R matrix
    evs <- eigen(R)$values
    abs.evs <- abs(evs)
    
  }
  
  k <- length(evs)
  
  if (method == "nyholt") {
    
    # effective number of tests (based on Nyholt, 2004)
    m <- 1 + (k - 1) * (1 - var(evs) / k)
    
  }
  
  if (method == "liji") {
    
    # effective number of tests (based on Li & Ji, 2005)
    # adding a small value to the eigenvalues to overcome numerical imprecisions
    abs.evs <- abs.evs + sqrt(.Machine$double.eps)
    m <- sum(ifelse(abs.evs >= 1, 1, 0) + (abs.evs - floor(abs.evs)))
    
  }
  
  if (method == "gao") {
    
    # effective number of tests (based on Gao, 2008)
    m <- which(cumsum(sort(abs.evs, decreasing = TRUE)) / sum(abs.evs) > 0.995)[1]
    
  }
  
  if (method == "galwey") {
    
    # effective number of tests (based on Galwey, 2009)
    evs[evs < 0] <- 0
    m <- sum(sqrt(evs))^2 / sum(evs)
    
  }
  
  # always round down estimated value
  m <- floor(m)
  
  return(m)
  
}
