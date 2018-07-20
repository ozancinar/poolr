fisher <- function(p, adjust = "none", m, R, size = 10000, seed, type = 2, 
                   emp.loop = FALSE, ...) {
  
  k <- length(p)
  
  # if m is provided by the user, then we don't need to check the adjustment method.
  if(!missing(m)) {
    
    m <- m
    adjust <- paste0(m, " (user defined)")
    
    testStat <- -2 * sum(log(p))
    testStat <- testStat * (m / k)
    pooled.p <- pchisq(testStat, df = 2 * m, lower.tail = FALSE)
    
    # warning the user if the user-defined m is larger than the number of p-values.
    if(m > k)
      warning("the user-defined effective number of test is larger than the number of p-values that were combined.")
    
  } else {
    
    # first, if the adjust is not given, it will be set to "none".
    if(missing(adjust))
      adjust <- "none"
    
    # now, checking the adjust argument.
    if(!adjust %in% c("none", "nyholt", "liji", "gao", "galwey", "empirical", "brown"))
      stop("adjust argument is not given correctly. Please refer to ?fisher for the correct set for adjust arguments.")
      
    if (adjust == "none") {
      
      testStat <- -2 * sum(log(p))
      pooled.p <- pchisq(testStat, df = 2*k, lower.tail = FALSE)
      adjust <- "none"
      
    } else if (adjust %in% c("nyholt", "liji", "gao", "galwey")) {
      
      m <- meff(R = R, method = adjust)
      adjust <- paste0("meff (", adjust, ")")
      
      testStat <- -2 * sum(log(p))
      testStat <- testStat * (m / k)
      pooled.p <- pchisq(testStat, df = 2 * m, lower.tail = FALSE)
      
    } else if (adjust == "brown") {
      
      tmp <- list(...)
      
      # if the Brown correlations are not provided by the user, we will use brown()
      # to get the Brown correlations.
      if (is.null(tmp$brownCov)) {
        
        covs <- brown(R)
        
      } else { # otherwise, the function will use user-defined correlations.
        
        covs <- tmp$brownCov
        
      }
      
      expx2 <- 2 * k
      varx2 <- sum(covs)
      fval <- 2 * expx2^2 / varx2
      cval <- varx2 / (2 * expx2)
      
      testStat <- -2 * sum(log(p))
      testStat <- testStat/cval
      pooled.p <- pchisq(testStat, df = fval, lower.tail = FALSE)
      adjust <- "brown"
      
    } else if (adjust == "empirical") {
      
      tmp <- list(...)
      
      # if an empirical distribution is not provided by the user, we will use 
      # empirical() to generate an empirical distribution.
      if (is.null(tmp$emp.dis)) {
        
        emp.dist <- empirical(R = R, method = "fisher", type = type, size = size, 
                              seed = seed, emp.loop = emp.loop)
        
      } else { # otherwise, the function will use the user-given empirical distribution.
        
        emp.dist <- tmp$emp.dist
        
      }
      
      testStat <- -2 * sum(log(p))
      pooled.p <- (sum(emp.dist >= testStat) + 1) / (size + 1)
      adjust <- "empirical"
      
    }
    
  }
  
  res <- list(p = pooled.p, testStat = testStat, adjust = adjust)
  class(res) <- "combP"
  return(res)

}
