
fisher <- function(p, adjust = "none", pca.method = NULL, R = NULL, size = 10000, seed = NULL, ...) {
   if(adjust == "none") {
      k <- length(p)
      testStat <- -2*sum(log(p))
      pooled.p <- pchisq(-2*sum(log(p)), df=2*k, lower.tail=FALSE)
   } else if(adjust == "m.eff") {
      k <- length(p)
      if(is.numeric(pca.method)) {
         eff <- pca.method 
      } else {
         eff <- meff(R = R, method = pca.method)
      }
      testStat <- -2 * sum(log(p)) * (eff / k)
      pooled.p <- pchisq(-2 * sum(log(p)) * (eff / k), df = 2 * eff, lower.tail = FALSE)
   } else if (adjust == "brown") {
      k <- length(p)
      
      tmp <- list(...)
      
      if(is.null(tmp$brownCov)) {
         covs <- brown(R)
      } else {
         covs <- tmp$brownCov
         if(is.vector(covs)) {
            if(length(covs) != k * (k - 1) / 2) stop("Length of the covariance vector is not correct.")
         } else if(is.matrix(covs)) {
            if(dim(covs) != k) stop("Dimensions of the covariance matrix is not correct.")
            covs <- covs[lower.tri(covs)]
         }
      }
      
      chi2val <- -2 * sum(log(p))
      expx2 <- 2 * k
      varx2 <- 4 * k + 2 * sum(covs)
      fval <- 2 * expx2^2 / varx2
      cval <- varx2 / (2 * expx2)
      
      testStat <- chi2val/cval
      pooled.p <- pchisq(chi2val/cval, df=fval, lower.tail=FALSE)
   } else if (adjust == "empirical") {
      k <- length(p)
      tmp.p <- pchisq(-2*sum(log(p)), df=2*k, lower.tail=FALSE)
      testStat <- -2*sum(log(p))
      method <- "fisher"
      
      tmp <- list(...)
      if(is.null(tmp$emp.dis)) {
         emp.dist <- empirical(p = p, R = R, method = method, size = size, seed = seed)
      } else {
         emp.dist <- tmp$emp.dist
      }
      
      pooled.p <- sum(emp.dist <= testStat) / length(emp.dist)
   }
   
   res <- list(p = pooled.p, testStat = testStat, adjust = paste0(adjust, " ", pca.method))
   class(res) <- "combP"
   return(res)
}
