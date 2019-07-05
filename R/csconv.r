csconv <- function(p, target) {
  
  test_fun <- NULL
  test_fun <- try(sys.call(-1))
  
  if(is.null(test_fun)) {
    if(missing(target)) {
      stop("target is missing")
    } else {
      target <- target
    }
  } else {
    if(missing(target)) {
      fun_info <- as.character(test_fun)
      fun_name <- fun_info[which(fun_info %in% c("fisher", "stouffer", "invchisq"))]
      
      if(fun_name == "fisher") {
        target <- "m2lp"
      } else if (fun_name == "stouffer") {
        target <- "z"
      } else if (fun_name == "invchisq") {
        target <- "chisq1"
      } else {
        stop("mvnconv() was called from a function that is not provided by poolR. This is not available for now.")
      }
    } else {
      target <- target
    }
  }
  
  if (target == "m2lp") {
    s <- -2 * log(p)
    rho_est <- max(0, 1 - var(s) / 4)
    k <- length(p)
    covs <- matrix(rho_est * 4, k, k)
    diag(covs) <- 4
  } else if (target == "z") {
    s <- qnorm(p)
    rho_est <- max(0, 1 - var(s) / 1)
    k <- length(p)
    covs <- matrix(rho_est * 1, k, k)
    diag(covs) <- 1
  } else if (target == "chisq1") {
    s <- qchisq(p, df = 1)
    rho_est <- max(0, 1 - var(s) / 2)
    k <- length(p)
    covs <- matrix(rho_est * 2, k, k)
    diag(covs) <- 2
  }
  
  return(covs)
}
