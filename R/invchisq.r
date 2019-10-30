invchisq <- function(p, adjust = "none", m, R, size = 10000, seed, side = 2,
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

   # set some defaults
   ci <- NULL
   m <- k

   if (adjust == "none") {

      pval <- pchisq(statistic, df = k, lower.tail = FALSE)
      attr(statistic, "df") <- k

   }

   if (adjust %in% c("nyholt", "liji", "gao", "galwey", "user")) {

      if (adjust != "user") {
         m <- meff(R = R, method = adjust)
      } else {
         # warn the user if the user-defined m is larger than the number of p-values
         if (m > k)
            warning("User-defined effective number of tests is larger than the number of p-values.")
      }

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

      tmp <- list(...)

      # checks/fixes for 'emp.step' argument
      emp.step <- .check.emp.step(emp.step, size = size, tmp = tmp)

      # observed pooled p-value
      pval.obs <- pchisq(statistic, df = k, lower.tail = FALSE)

      for (i in 1:length(emp.step$size)) {

         if (!is.null(tmp$verbose) && tmp$verbose)
            cat("Size:", emp.step$size[i], " Threshold:", emp.step$thres[i], "\n")

         size <- emp.step$size[i]

         if (is.null(tmp$emp.dist)) {
            emp.dist <- empirical(R = R, method = fun, side = side,
                                  size = size, seed = seed, emp.loop = emp.loop)
         } else {
            emp.dist <- tmp$emp.dist
         }

         pval <- (sum(emp.dist <= pval.obs) + 1) / (size + 1)

         if (pval >= emp.step$thres[i]) {
            attr(statistic, "df") <- k
            ci <- as.numeric(binom.test((sum(emp.dist <= pval.obs) + 1), (size + 1))$conf.int)
            break
         }

      }

   }

   res <- list(p = pval, ci = ci, k = k, m = m, fun = fun, adjust = adjust, statistic = statistic)

   class(res) <- "combp"
   return(res)

}
