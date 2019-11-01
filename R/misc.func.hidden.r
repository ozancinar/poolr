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

}

.check.m <- function(R, adjust, m, k) {

   if (adjust != "user") {
      m <- meff(R = R, method = adjust)
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

   # check if 'size' is numeric
   if (!is.numeric(size)) 
     stop("Argument 'size' must be numeric. See help(", call.fun, ").", call.=FALSE)

   # check if 'size' contains appropriate values
   if (any(size < 1)) {
     stop("Argument 'size' must include values greater than 1. See help(", call.fun, ").", call.=FALSE) 
   }

   # check if 'thres' is appropriate when it is given
   if (!is.null(thres)) {
      if (!is.numeric(thres))
      stop("Argument 'thres' must be numeric. See help(", call.fun, ").", call.=FALSE)

      if (any(thres > 1) || any(thres < 0))
     stop("Argument 'thres' must include values between 0 and 1. See help(", call.fun, ").", call.=FALSE) 
   }

   # check if 'thres' is given when 'size' is a vector
   if (length(size) > 1 & is.null(thres)) {
      stop("Argument 'thres' must be specified when 'size' is a vector. See help(", call.fun, ").", call.=FALSE)
   }

   # check the compatibility of the length of 'thres' when 'size' is a vector
   if (length(size) > 1) {
      step.length <- length(size)
      compatible.lengths <- c(1, step.length, step.length - 1)
      if (!length(thres) %in% compatible.lengths) 
         stop("Incompatible vector length for 'thres'. See help(", call.fun, ").", call.=FALSE)
   }

   # if 'size' is set to a vector while an empirical distribution is provided by the user
   if (length(size) > 1 & !is.null(ddd$emp.dist))
      stop("Stepwise algorithm cannot be used with a user-defined empirical distribution. See help(", call.fun, ").", call.=FALSE)

   # if 'thres' is set to a vector while 'size' is a single numeric
   if (length(size) == 1 & length(thres) > 1) {
      emp.setup <- list(size = size, thres = 0)
      warning("Multiple thresholds cannot be used with a single sample size. Threshold was set to 0. See help(", call.fun, ").", call.=FALSE)
   }

   # if 'size' is a single numeric, set it to the chosen size with threshold 0
   if (length(size) == 1 || !is.null(ddd$emp.dist))
      emp.setup <- list(size = size, thres = 0)

   # if 'thres' is a single numeric when 'size' is a vector
   if (length(size) > 1 & length(thres) == 1) {
      emp.setup <- list(size = size, thres = c(rep(thres, length(size) - 1), 0))
   }

   # if the length of 'thres' is one less than the length of 'size'
   if (length(size) > 1 & length(thres) == length(size) - 1) {
      emp.setup <- list(size = size, thres = c(thres, 0))
   }

   # if 'thres' and 'size' have the same length
   if (length(size) > 1 & length(thres) == length(size)) {
      emp.setup <- list(size = size, thres = c(thres[1:length(size) - 1], 0))
      warning("Arguments 'thres' and 'size' have the same length. The last threshold is ignored. See help(", call.fun, ").", call.=FALSE)
   }

   return(emp.setup)

}

.do.emp <- function(pval.obs, emp.setup, ddd, R, method, side, emp.loop) {

   for (i in 1:length(emp.setup$size)) {

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
