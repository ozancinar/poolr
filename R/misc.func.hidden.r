############################################################################

.check.p <- function(p) {

   if (missing(p))
      stop("Argument 'p' must be specified.", call.=FALSE)

   if (!is.vector(p) || !is.numeric(p))
      stop("Argument 'p' must be a numeric vector.", call.=FALSE)

   if (any(is.na(p)))
      stop("Values in 'p' vector must not contain NAs.", call.=FALSE)

   if (any(p < 0) || any(p > 1))
      stop("Values in 'p' vector (i.e., the p-values) must be between 0 and 1.", call.=FALSE)

}

.check.R <- function(R, checksym=TRUE, checkna=TRUE, checkpd=FALSE, checkcor=FALSE, checkdiag = 1, isbase=TRUE, k, adjust, fun) {

   # turn a "dpoMatrix" object (from nearPD()) into a 'plain' matrix (since is.matrix() is FALSE for such objects)
   if (inherits(R, "dpoMatrix"))
      R <- as.matrix(R)

   # if R is a dataframe, turn it into a matrix
   if (inherits(R, "data.frame"))
      R <- as.matrix(R)

   # check that 'R' is a symmetric matrix
   if (checksym && !(is.matrix(R) && isSymmetric(unname(R))))
      stop("Argument 'R' must be a symmetric matrix.", call.=FALSE)

   # check if 'R' contains NAs
   if (checkna && any(is.na(R)))
      stop("Values in 'R' vector must not contain NAs.", call.=FALSE)

   # check if 'R' is positive definite; if not, make it
   if (checkpd && any(eigen(R)$values <= 0)) {
      R <- as.matrix(Matrix::nearPD(R, corr=TRUE)$mat)
      warning("Matrix 'R' is not positive definite. Used Matrix::nearPD() to make 'R' positive definite.", call.=FALSE)
   }

   # check that all values in R are between -1 and 1
   if (checkcor && any(abs(R) > 1))
      stop("Argument 'R' must be a correlation matrix, but contains values outside [-1,1].", call.=FALSE)

   # check that all diagonal values are equal to 1
   if (length(unique(diag(R))) != 1 || unique(diag(R)) != checkdiag)
      stop("Argument 'R' must be a correlation matrix with diagonal values are all equal to 1.", call.=FALSE)

   # checks that are relevant only when called from the base functions

   if (isbase) {

      # check that dimensions of 'R' match the length of 'p'
      if (k != nrow(R))
         stop("Length of 'p' vector (", k, ") does not match the dimensions of the 'R' matrix (", nrow(R), ",", ncol(R), ").", call.=FALSE)

      # check if user specified 'R' argument but no adjustment method
      if (adjust == "none")
         warning("Although argument 'R' was specified, no adjustment method was chosen via the 'adjust' argument.\n  To account for dependence, specify an adjustment method. See help(", fun, ") for details.", call.=FALSE)

      # if 'm' has been specified, then warn the user that 'R' matrix is actually ignored
      if (adjust == "user")
         warning("When 'm' is specified, argument 'R' is irrelevant and ignored.")

   }

   return(R)

}

.check.m <- function(R, adjust, m, k, ...) {

   if (adjust != "user") {
      m <- meff(R = R, method = adjust, ...)
   } else {
      # warn the user if the user-defined m is larger than the number of p-values
      if (m > k)
         warning("User-defined effective number of tests is larger than the number of p-values.")
   }

   return(m)

}

.check.side <- function(side) {

   if (length(side) != 1)
      stop("Argument 'side' must be of length 1.")

   if (!(side %in% c(1,2)))
      stop("Argument 'side' must be either 1 or 2.")

}

############################################################################

.check.emp.setup <- function(size, threshold, ddd) {

   # get name of calling function
   call.fun <- as.character(sys.call(-1)[1])

   if (!is.null(ddd$emp.dist)) {

      # special case: empirical distribution is provided by the user

      # if 'size' is set to a vector or 'threshold' is not null, issue a warning
      if (length(size) > 1L || !is.null(threshold))
         warning("Stepwise algorithm cannot be used with a user-defined empirical distribution. See help(", call.fun, ").", call.=FALSE)

      # set 'size' to the length of the empirical distribution provided by the user
      size <- length(ddd$emp.dist)

      # and set 'threshold' to 0
      emp.setup <- list(size = size, threshold = 0)

      return(emp.setup)

   }

   # check if 'size' is numeric
   if (!is.numeric(size))
     stop("Argument 'size' must be numeric. See help(", call.fun, ").", call.=FALSE)

   # check if all values in 'size' are >= 1
   if (any(size < 1))
     stop("Values in 'size' must be >= 1. See help(", call.fun, ").", call.=FALSE)

   #                                                              threshold
   #             | NULL                        | scalar                                          | vector
   #      -------+-----------------------------+-------------------------------------------------+------------------------------------------
   #      scalar | ok, set threshold=0         | set threshold=0 and warn that this happened     | <- same
   # size -------+-----------------------------+-------------------------------------------------+------------------------------------------
   #      vector | error (must spec threshold) | ok, rep threshold as needed, set last to 0      | ok, check lengths and make sure last is 0
   #      -------+-----------------------------+-------------------------------------------------+------------------------------------------

   length.size  <- length(size)

   if (length.size == 1L) {

      # if 'size' is a scalar, argument 'threshold' is ignored (it is always set to 0)
      if (!is.null(threshold))
         warning("Argument 'threshold' not relevant with a single sample size. Threshold was set to 0. See help(", call.fun, ").", call.=FALSE)

      threshold <- 0

   } else {

      # if 'size' is a vector but 'threshold' was not specified, issue an error
      if (is.null(threshold))
         stop("Argument 'threshold' must be specified when 'size' is a vector. See help(", call.fun, ").", call.=FALSE)

      # check if 'threshold' is appropriate when it is given
      if (!is.numeric(threshold))
         stop("Argument 'threshold' must be numeric. See help(", call.fun, ").", call.=FALSE)
      if (any(threshold > 1) || any(threshold < 0))
         stop("Values in 'threshold' must be between 0 and 1. See help(", call.fun, ").", call.=FALSE)

      # if 'threshold' is a scalar, rep it as needed and set last value to 0
      if (length(threshold) == 1L)
         threshold <- c(rep(threshold, length.size - 1), 0)

      # if length of 'threshold' is one less than length of 'size', add 0 to the 'threshold' vector
      if (length(threshold) == length.size - 1)
         threshold <- c(threshold, 0)

      # check compatibility of the lengths of the 'size' and 'threshold' arguments
      if (length(threshold) != length.size)
         stop("Length of 'threshold' argument is not compatible with length of 'size' argument. See help(", call.fun, ").", call.=FALSE)

      length.threshold <- length(threshold)

      if (threshold[length.threshold] != 0) {
         warning("Last value of 'threshold' is not 0. Last threshold set to 0. See help(", call.fun, ").", call.=FALSE)
         threshold[length.threshold] <- 0
      }

   }

   emp.setup <- list(size = size, threshold = threshold)

   return(emp.setup)

}

.do.emp <- function(pval.obs, emp.setup, ddd, R, method, side, batchsize) {

   for (j in seq_along(emp.setup$size)) {

      if (!is.null(ddd$verbose) && ddd$verbose)
         cat("Size:", emp.setup$size[j], " Threshold:", emp.setup$threshold[j], "\n")

      size <- emp.setup$size[j]

      if (is.null(ddd$emp.dist)) {
         empirical.args <- c(list(R = R, method = method, side = side, size = size, batchsize = batchsize), ddd)
         emp.dist <- do.call(empirical, empirical.args)
         #emp.dist <- empirical(R = R, method = method, side = side, size = size, batchsize = batchsize)
      } else {
         emp.dist <- ddd$emp.dist
      }

      pval <- (sum(emp.dist <= pval.obs) + 1) / (size + 1)

      if (pval >= emp.setup$threshold[j]) {
         ci <- as.numeric(binom.test((sum(emp.dist <= pval.obs) + 1), (size + 1))$conf.int)
         break
      }

   }

   return(list(pval = pval, ci = ci))

}

############################################################################

# simplified version of MASS::mvrnorm()

.simmvn <- function(n = 1, Sigma, method = "mvt_eigen") {

   p <- nrow(Sigma)
   eS <- eigen(Sigma, symmetric = TRUE)
   eval <- eS$values
   evec <- eS$vectors

   X <- matrix(rnorm(p * n), nrow = n, byrow = TRUE)

   if (method == "mass_eigen") {   
      return(X %*% diag(sqrt(pmax(eval, 0)), p) %*% t(evec))
   } else if (method == "mvt_eigen") {
      return(X %*% t(evec %*% (t(evec) * sqrt(pmax(eval, 0)))))
   }

}

############################################################################

.fisher <- function(p, k, alpha) {
   statistic <- -2 * sum(log(p))
   pchisq(statistic, df = 2 * k, lower.tail = FALSE)
}

.stouffer <- function(p, k, alpha) {
   statistic <- sum(qnorm(p, lower.tail = FALSE)) / sqrt(k)
   pnorm(statistic, lower.tail = FALSE)
}

.invchisq <- function(p, k, alpha) {
   statistic <- sum(qchisq(p, df = 1, lower.tail = FALSE))
   pchisq(statistic, df = k, lower.tail = FALSE)
}

.bonferroni <- function(p, k, alpha) {
   statistic <- min(p)
   min(1, statistic * k)
}

.tippett <- function(p, k, alpha) {
   statistic <- min(p)
   1 - (1 - statistic)^k
}

.binotest <- function(p, k, alpha) {
   statistic <- sum(p <= alpha)
   sum(dbinom(statistic:k, k, alpha))
}

############################################################################
