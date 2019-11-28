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

.check.R <- function(R, k, adjust, fun) {

   # check that 'R' is a symmetric matrix
   if (!is.matrix(R) || !isSymmetric(unname(R)))
      stop("Argument 'R' must be a symmetric matrix.", call.=FALSE)

   # check if 'R' contains NAs
   if (any(is.na(R)))
      stop("Values in 'R' vector must not contain NAs.", call.=FALSE)

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

.check.emp.setup <- function(size, thres, ddd) {

   # get name of calling function
   call.fun <- as.character(sys.call(-1)[1])

   if (!is.null(ddd$emp.dist)) {

      # special case: empirical distribution is provided by the user

      # if 'size' is set to a vector or 'thres' is not null, issue a warning
      if (length(size) > 1L || !is.null(thres))
         warning("Stepwise algorithm cannot be used with a user-defined empirical distribution. See help(", call.fun, ").", call.=FALSE)

      # set 'size' to the length of the empirical distribution provided by the user
      size <- length(ddd$emp.dist)

      # and set 'thres' to 0
      emp.setup <- list(size = size, thres = 0)

      return(emp.setup)

   }

   # check if 'size' is numeric
   if (!is.numeric(size))
     stop("Argument 'size' must be numeric. See help(", call.fun, ").", call.=FALSE)

   # check if all values in 'size' are >= 1
   if (any(size < 1))
     stop("Values in 'size' must be >= 1. See help(", call.fun, ").", call.=FALSE)

   #                                     thres
   #             | NULL                    | scalar                                  | vector
   #      -------+-------------------------+-----------------------------------------+------------------------------------------
   #      scalar | ok, set thres=0         | set thres=0 and warn that this happened | <- same
   # size -------+-------------------------+-----------------------------------------+------------------------------------------
   #      vector | error (must spec thres) | ok, rep thres as needed, set last to 0  | ok, check lengths and make sure last is 0
   #      -------+-------------------------+-----------------------------------------+------------------------------------------

   length.size  <- length(size)

   if (length.size == 1L) {

      # if 'size' is a scalar, argument 'thres' is ignored (it is always set to 0)
      if (!is.null(thres))
         warning("Argument 'thres' not relevant with a single sample size. Threshold was set to 0. See help(", call.fun, ").", call.=FALSE)

      thres <- 0

   } else {

      # if 'size' is a vector but 'thres' was not specified, issue an error
      if (is.null(thres))
         stop("Argument 'thres' must be specified when 'size' is a vector. See help(", call.fun, ").", call.=FALSE)

      # check if 'thres' is appropriate when it is given
      if (!is.numeric(thres))
         stop("Argument 'thres' must be numeric. See help(", call.fun, ").", call.=FALSE)
      if (any(thres > 1) || any(thres < 0))
         stop("Values in 'thres' must be between 0 and 1. See help(", call.fun, ").", call.=FALSE)

      # if 'thres' is a scalar, rep it as needed and set last value to 0
      if (length(thres) == 1L)
         thres <- c(rep(thres, length.size - 1), 0)

      # if length of 'thres' is one less than length of 'size', add 0 to the 'thres' vector
      if (length(thres) == length.size - 1)
         thres <- c(thres, 0)

      # check compatibility of the lengths of the 'size' and 'thres' arguments
      if (length(thres) != length.size)
         stop("Length of 'thres' argument is not compatible with length of 'size' argument. See help(", call.fun, ").", call.=FALSE)

      length.thres <- length(thres)

      if (thres[length.thres] != 0) {
         warning("Last value of 'thres' is not 0. Last threshold set to 0. See help(", call.fun, ").", call.=FALSE)
         thres[length.thres] <- 0
      }

   }

   emp.setup <- list(size = size, thres = thres)

   return(emp.setup)

}

.do.emp <- function(pval.obs, emp.setup, ddd, R, method, side, emp.loop) {

   for (i in seq_along(emp.setup$size)) {

      if (!is.null(ddd$verbose) && ddd$verbose)
         cat("Size:", emp.setup$size[i], " Threshold:", emp.setup$thres[i], "\n")

      size <- emp.setup$size[i]

      if (is.null(ddd$emp.dist)) {
         emp.dist <- empirical(R = R, method = method, side = side,
                               size = size, emp.loop = emp.loop)
      } else {
         emp.dist <- ddd$emp.dist
      }

      pval <- (sum(emp.dist <= pval.obs) + 1) / (size + 1)

      if (pval >= emp.setup$thres[i]) {
         ci <- as.numeric(binom.test((sum(emp.dist <= pval.obs) + 1), (size + 1))$conf.int)
         break
      }

   }

   return(list(pval = pval, ci = ci))

}