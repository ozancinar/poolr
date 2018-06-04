meff <- function(x, eigen = FALSE, method = NULL, ...) {

   # Check if the method is correct.
   if(!method %in% c("nyholt", "li.ji", "gao", "galwey")) {
         stop("The method for PCA is not correct.")
   }

   # if the input x is eigen values, i.e., eigen is TRUE, then it will not calculate the eigen values of the matrix.
   if(eigen) {
      if(!class(x) %in% c("numeric", "integer")) {
         stop("eigen values are not numeric or integer.")
      }
      k <- length(x)

      evs <- x
      abs.evs <- abs(evs)
   } else {
      # number of p-values
      k <- nrow(x)

      # dimension checks
      if(!isSymmetric(x)) {
         stop("R is not symmetric.")
      }

      # Converting the correlation matrix to positive-definite.
      x <- as.matrix(nearPD(x)$mat)

      ### get eigenvalues and absolute eigenvalues of R matrix
      evs <- eigen(x)$values
      abs.evs <- abs(evs)
   }


   if(method == "nyholt") {
      ### effective number of tests (based on Nyholt, 2004)
      eff <- 1 + (k - 1) * (1 - var(evs) / k)
      eff <- floor(eff)
   } else if(method == "li.ji") {
      ### effective number of tests (based on Li & Ji, 2005)
      # adding a small value to the eigen values in order to overcome the numeric calculation problem.
      abs.evs <- abs.evs + sqrt(.Machine$double.eps)
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
