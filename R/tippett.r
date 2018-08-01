tippett <- function(p, adjust = "none", m, R, size = 10000, seed, type = 2, 
                    emp.loop = FALSE, ...) {
  
  k <- length(p)
  
  # if m is provided by the user, then we don't need to check the adjustment method.
  if(!missing(m)) {
    
    m <- m
    adjust <- paste0(m, " (user defined)")
    
    testStat <- 1 - (1 - min(p))^m
    pooled.p <- testStat
    
    # warning the user if the user-defined m is larger than the number of p-values.
    if(m > k)
      warning("the user-defined effective number of test is larger than the number of p-values that were combined.")
    
  } else {
    
    # first, if the adjust is not given, it will be set to "none".
    if(missing(adjust))
      adjust <- "none"
    
    # now, checking the adjust argument.
    if(!adjust %in% c("none", "nyholt", "liji", "gao", "galwey", "empirical"))
      stop("adjust argument is not given correctly. Please see ?tippett for the correct options for adjust.")
      
    if (adjust == "none") {
      
      testStat <- 1 - (1 - min(p))^k
      pooled.p <- testStat
      adjust <- "none"
      
    } else if (adjust %in% c("nyholt", "liji", "gao", "galwey")) {
      
      m <- meff(R = R, method = adjust)
      adjust <- paste0("meff (", adjust, ")")
      
      testStat <- 1 - (1 - min(p))^m
      pooled.p <- testStat
      
    } else if (adjust == "empirical") {
      
      tmp <- list(...)
      
      # if an empirical distribution is not provided by the user, we will use 
      # empirical() to generate an empirical distribution.
      if (is.null(tmp$emp.dis)) {
        
        emp.dist <- empirical(R = R, method = "tippett", type = type, size = size, 
                              seed = seed, emp.loop = emp.loop)
        
      } else { # otherwise, the function will use the user-given empirical distribution.
        
        emp.dist <- tmp$emp.dist
        
      }
      
      testStat <- 1 - (1 - min(p))^k
      pooled.p <- (sum(emp.dist <= testStat) + 1) / (size + 1)
      adjust <- "empirical"
      
    }
    
  }
  
  res <- list(p = pooled.p, testStat = testStat, adjust = adjust)
  class(res) <- "combP"
  return(res)

}
