invchisq <- function(p, adjust = "none", m, R, size = 10000, side = 2,
                     emp.loop = FALSE, emp.step, ...) {

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
   statistic <- sum(qchisq(p, df = 1, lower.tail = FALSE))
   attr(statistic, "df") <- k

   # set some defaults
   ci <- NULL
   m <- k

   if (adjust == "none") {

      pval <- pchisq(statistic, df = k, lower.tail = FALSE)

   }

   if (adjust %in% c("nyholt", "liji", "gao", "galwey", "user")) {

      m <- .check.m(R = R, adjust = adjust, m = m, k = k)

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

      ddd <- list(...)

      # checks/fixes for 'emp.step' argument
      emp.step <- .check.emp.step(emp.step, size = size, ddd = ddd)

      # observed pooled p-value
      pval.obs <- pchisq(statistic, df = k, lower.tail = FALSE)

      # get empirically derived p-value
      tmp <- .do.emp(pval.obs = pval.obs, emp.step = emp.step, ddd = ddd,
                     R = R, method = fun, side = side, emp.loop = emp.loop)

      pval <- tmp$pval
      ci <- tmp$ci

   }

   res <- list(p = pval, ci = ci, k = k, m = m, fun = fun, adjust = adjust, statistic = statistic)

   class(res) <- "combp"
   return(res)

}
