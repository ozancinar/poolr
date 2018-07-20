binotest <- function(p, adjust = "none", m, R, alpha = 0.05, size = 10000, seed, 
                     type = 2, ...) {
  
  k <- length(p)
  r <- sum(p <= alpha)
  
  # if m is provided by the user, then we don't need to check the adjustment method.
  if(!missing(m)) {
    
    m <- m
    adjust <- paste0(m, " (user defined)")
    
    testStat <- dbinom(round(r * m / k), m, alpha)
    pooled.p <- sum(dbinom(round(r * m / k):m, m, alpha))
    
  } else {
    
    # if m is not provided by the user, then the functions will use adjust argument.
    # checking if the adjust argument is given correctly.
    
    # first, if the adjust is not given, it will be set to "none".
    if(missing(adjust)) {
      
      adjust <- "none"
      
    }
    
    # now, checking the adjust argument.
    
    if(!adjust %in% c("none", "nyholt", "liji", "gao", "galwey", "empirical")) {
      
      stop("adjust argument is not given correctly. Please refer to ?binotest for the correct set for adjust arguments.")
      
    }
    
    if (adjust == "none") {
      
      testStat <- dbinom(r, k, alpha)
      pooled.p <- sum(dbinom(r:k, k, alpha))
      adjust <- "none"
      
    } else if (adjust %in% c("nyholt", "liji", "gao", "galwey")) {
      
      m <- meff(R = R, method = adjust)
      adjust <- paste0("meff (", adjust, ")")
      
      testStat <- dbinom(round(r * m / k), m, alpha)
      pooled.p <- sum(dbinom(round(r * m / k):m, m, alpha))
      
    } else if (adjust == "empirical") {
      
      tmp <- list(...)
      
      # if an empirical distribution is not provided by the user, we will use 
      # empirical() to generate an empirical distribution.
      if (is.null(tmp$emp.dis)) {
        
        emp.dist <- empirical(R = R, method = "binotest", type = type, size = size, 
                              seed = seed)
        
      } else { # otherwise, the function will use the user-given empirical distribution.
        
        emp.dist <- tmp$emp.dist
        
      }
      
      testStat <- dbinom(r, k, alpha)
      pooled.p <- (sum(emp.dist <= testStat) + 1) / (size + 1)
      adjust <- "empirical"
      
    }
    
  }
  
  res <- list(p = pooled.p, testStat = testStat, adjust = adjust)
  class(res) <- "combP"
  return(res)

}
