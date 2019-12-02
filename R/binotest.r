binotest <- function(p, adjust = "none", m, R, alpha = 0.05, size = 10000, thres, side = 2,
                     emploop = FALSE, ...) {

   # checks for 'p' argument
   .check.p(p)

   k <- length(p)

   # match 'adjust' argument
   adjust <- match.arg(adjust, c("none", "nyholt", "liji", "gao", "galwey", "empirical"))

   # if m is specified, apply effective number of test adjustment with user-defined m
   if (!missing(m))
      adjust <- "user"

   # get name of function
   fun <- as.character(sys.call()[1])

   if (missing(R)) {

      # check if 'R' is specified when using an adjustment method (does not apply to "user")
      if (adjust %in% c("nyholt", "liji", "gao", "galwey", "empirical", "generalized"))
         stop("Argument 'R' must be specified when using an adjustment method.")

   } else {

      R # force evaluation of 'R' argument, so that R=mvnconv(R) works

      # checks for 'R' argument
      .check.R(R, k = k, adjust = adjust, fun = fun)

   }

   # compute test statistic
   statistic <- sum(p <= alpha)

   # set some defaults
   ci <- NULL
   if (adjust != "user")
      m <- NULL

   if (adjust == "none") {

      pval <- sum(dbinom(statistic:k, k, alpha))

   }

   if (adjust %in% c("nyholt", "liji", "gao", "galwey", "user")) {

      m <- .check.m(R = R, adjust = adjust, m = m, k = k, ...)

      pval <- sum(dbinom(round(statistic * m / k):m, m, alpha))

   }

   if (adjust == "empirical") {

      ddd <- list(...)

      # setting 'thres' to NULL if it is missing for further checks
      if (missing(thres))
         thres <- NULL

      # checks/fixes for 'size' and 'thres' arguments
      emp.setup <- .check.emp.setup(size = size, thres = thres, ddd = ddd)

      # observed pooled p-value
      pval.obs <- sum(dbinom(statistic:k, k, alpha))

      # get empirically derived p-value
      tmp <- .do.emp(pval.obs = pval.obs, emp.setup = emp.setup, ddd = ddd,
                     R = R, method = fun, side = side, emploop = emploop)

      pval <- tmp$pval
      ci <- tmp$ci

   }

   res <- list(p = pval, ci = ci, k = k, m = m, fun = fun, adjust = adjust, statistic = statistic)

   class(res) <- "combp"
   return(res)

}
