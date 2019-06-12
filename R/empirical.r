
empirical <- function(R, method, type, size = 10000, seed, emp.loop, ...) {
  
  # match method argument
  method <- match.arg(method, c("fisher", "stouffer", "invchisq", "binotest", 
                                "bonferroni", "tippett"))

  # check type argument
  if(!type %in% c(1, 2))
      stop("argument 'type' must be either 1 or 2.")
  
  if(missing(emp.loop))
    emp.loop <- FALSE

  k <- nrow(R)

  # check that R is symmetric
  if(!isSymmetric(R))
      stop("R is not symmetric.")

  # checking if the correlation matrix is positive definite by testing if all 
  # eigen-values are positive.
  if(!all(eigen(R)$values > 0)) {
    R <- nearPD(R)$mat
    warning("the correlation matrix is not positive definite. empirical() used Matrix::nearPD to make the correlation matrix positive definite.")
  }
  
  if(!missing(seed))
    set.seed(seed)
  
  if(emp.loop == FALSE) {
    z <- mvrnorm(size, mu = rep(0, k), Sigma = R)
    
    if(type == 1)
      p <- pnorm(z, lower.tail = FALSE)
    if(type == 2)
      p <- 2 * pnorm(abs(z), lower.tail = FALSE)
    
    if(method == "fisher")
      emp_dist <- apply(p, 1, function(x) fisher(x)$test_stat)
    
    if(method == "stouffer")
      emp_dist <- apply(p, 1, function(x) stouffer(x)$test_stat)
    
    if(method == "invchisq")
      emp_dist <- apply(p, 1, function(x) invchisq(x)$test_stat)
    
    if(method == "binotest")
      emp_dist <- apply(p, 1, function(x) binotest(x)$test_stat)
    
    if(method == "bonferroni")
      emp_dist <- apply(p, 1, function(x) bonferroni(x)$test_stat)
    
    if(method == "tippett")
      emp_dist <- apply(p, 1, function(x) tippett(x)$test_stat)
    
  } else {
    
    emp_dist <- NULL
    
    for(i in 1:size) {
      z <- mvrnorm(1, mu = rep(0, k), Sigma = R)
      
      if(type == 1)
        p <- pnorm(z, lower.tail = FALSE)
      if(type == 2)
        p <- 2 * pnorm(abs(z), lower.tail = FALSE)
      
      if(method == "fisher")
        emp_dist <- append(emp_dist, fisher(p)$test_stat)
      
      if(method == "stouffer")
        emp_dist <- append(emp_dist, stouffer(p)$test_stat)
      
      if(method == "invchisq")
        emp_dist <- append(emp_dist, invchisq(p)$test_stat)
      
      if(method == "binotest")
        emp_dist <- append(emp_dist, binotest(p)$test_stat)
      
      if(method == "bonferroni")
        emp_dist <- append(emp_dist, bonferroni(p)$test_stat)
      
      if(method == "tippett")
        emp_dist <- append(emp_dist, tippett(p)$test_stat)
    }
  }
  return(emp_dist)
}
