invchisq <- function(p, adjust = "none", R, m, size = 10000, threshold, side = 2, batchsize, ...) {

   # checks for 'p' argument
   .check.p(p)

   k <- length(p)

   # match 'adjust' argument
   adjust <- match.arg(adjust, c("none", "nyholt", "liji", "gao", "galwey", "empirical", "generalized"))

   # if m is specified, apply effective number of test adjustment with user-defined m
   if (!missing(m))
      adjust <- "user"

   # get name of function
   fun <- as.character(sys.call()[1])
   fun <- gsub("^poolr::", "", fun)

   ddd <- list(...)

   if (missing(R)) {

      # check if 'R' is specified when using an adjustment method (does not apply to "user")
      if (adjust %in% c("nyholt", "liji", "gao", "galwey", "empirical", "generalized"))
         stop("Argument 'R' must be specified when using an adjustment method.")

   } else {

      R # force evaluation of 'R' argument, so that R=mvnconv(R) works

      # checks for 'R' argument
      if (adjust == "generalized") {

         R <- .check.R(R, k = k, adjust = adjust, checkdiag = 2, fun = fun)
         
      } else {

         R <- .check.R(R, k = k, adjust = adjust, fun = fun)

      }

   }

   # compute test statistic
   statistic <- sum(qchisq(p, df = 1, lower.tail = FALSE))
   attr(statistic, "df") <- k

   # set some defaults
   ci <- NULL
   if (adjust != "user")
      m <- NULL

   if (adjust == "none") {

      pval <- pchisq(statistic, df = k, lower.tail = FALSE)

   }

   if (adjust %in% c("nyholt", "liji", "gao", "galwey", "user")) {

      m <- .check.m(R = R, adjust = adjust, m = m, k = k, ...)

      statistic <- statistic * (m / k)
      pval <- pchisq(statistic, df = m, lower.tail = FALSE)
      attr(statistic, "df") <- m

   }

   if (adjust == "generalized") {

      covs  <- R
      expx2 <- k
      varx2 <- sum(covs)
      fval  <- 2 * expx2^2 / varx2
      cval  <- varx2 / (2 * expx2)

      statistic <- statistic / cval
      pval <- pchisq(statistic, df = fval, lower.tail = FALSE)
      attr(statistic, "df") <- fval

   }

   if (adjust == "empirical") {

      # setting 'batchsize' to NULL if it is missing
      if (missing(batchsize))
         batchsize <- NULL

      # setting 'threshold' to NULL if it is missing for further checks
      if (missing(threshold))
         threshold <- NULL

      # checks/fixes for 'size' and 'threshold' arguments
      emp.setup <- .check.emp.setup(size = size, threshold = threshold, ddd = ddd)

      # observed pooled p-value
      pval.obs <- pchisq(statistic, df = k, lower.tail = FALSE)

      # get empirically derived p-value
      tmp <- .do.emp(pval.obs = pval.obs, emp.setup = emp.setup, ddd = ddd,
                     R = R, method = fun, side = side, batchsize = batchsize)

      pval <- tmp$pval
      ci <- tmp$ci

   }

   res <- list(p = c(pval), ci = ci, k = k, m = m, adjust = adjust, statistic = statistic, fun = fun)

   class(res) <- "poolr"
   return(res)

}
