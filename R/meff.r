meff <- function(R, eigen, method, ...) {

   # match 'method' argument
   method <- match.arg(method, c("nyholt", "liji", "gao", "galwey"))

   if (missing(eigen)) {

      # check if 'R' is specified
      if (missing(R))
         stop("Argument 'R' must be specified.", call.=FALSE)

      # check that 'R' is a symmetric matrix
      if (!is.matrix(R) || !isSymmetric(unname(R)))
         stop("Argument 'R' must be a symmetric matrix.", call.=FALSE)

      # ensure that the correlation matrix is positive semi-definite
      #R <- as.matrix(Matrix::nearPD(R)$mat)

      # get eigenvalues of 'R' matrix
      evs <- base::eigen(R)$values

   } else {

      # can pass eigenvalues directly to function via 'eigen'

      if (!is.vector(eigen) || !is.numeric(eigen))
         stop("Argument 'eigen' must be a numeric vector.", call.=FALSE)

      evs <- eigen

   }

   if (method == "nyholt") {

      # effective number of tests (based on Nyholt, 2004)
      k <- length(evs)
      m <- 1 + (k - 1) * (1 - var(evs) / k)

   }

   if (method == "liji") {

      # effective number of tests (based on Li & Ji, 2005)
      # adding a small value to the absolute eigenvalues to overcome numerical imprecisions
      abs.evs <- abs(evs) + sqrt(.Machine$double.eps)
      m <- sum(ifelse(abs.evs >= 1, 1, 0) + (abs.evs - floor(abs.evs)))

   }

   if (method == "gao") {

      # effective number of tests (based on Gao, 2008)

      ddd <- list(...)

      # allow user to specify value of C via ... but otherwise use 0.995
      if (!is.null(ddd$C)) {
         C <- ddd$C
      } else {
         C <- 0.995
      }

      m <- which(cumsum(sort(evs, decreasing = TRUE)) / sum(evs) > C)[1]

   }

   if (method == "galwey") {

      # effective number of tests (based on Galwey, 2009)
      evs[evs < 0] <- 0
      m <- sum(sqrt(evs))^2 / sum(evs)

   }

   # always round down the estimated value
   m <- floor(m)

   return(m)

}
