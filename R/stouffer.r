stouffer <- function(p, adjust = "none", m, R, size = 10000, seed, side = 2,
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
   statistic <- sum(qnorm(p, lower.tail = FALSE)) / sqrt(k)

   # set some defaults
   ci <- NULL
   m <- k

   if (adjust == "none") {

      pval <- pnorm(statistic, lower.tail = FALSE)

   }

   if (adjust %in% c("nyholt", "liji", "gao", "galwey", "user")) {

      if (adjust != "user") {
         m <- meff(R = R, method = adjust)
      } else {
         # warn the user if the user-defined m is larger than the number of p-values
         if (m > k)
            warning("User-defined effective number of tests is larger than the number of p-values.")
      }

      statistic <- statistic * sqrt(m / k)
      pval <- pnorm(statistic, lower.tail = FALSE)

   }

   if (adjust == "generalized") {

      statistic <- statistic * sqrt(k) / sqrt(sum(R))
      pval <- pnorm(statistic, lower.tail = FALSE)

   }

   if (adjust == "empirical") {

      tmp <- list(...)

      # checks/fixes for 'emp.step' argument
      emp.step <- .check.emp.step(emp.step, size = size, tmp = tmp)

      # observed pooled p-value
      pval.obs <- pnorm(statistic, lower.tail = FALSE)

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
            ci <- as.numeric(binom.test((sum(emp.dist <= pval.obs) + 1), (size + 1))$conf.int)
            break
         }

      }

   }

   res <- list(p = pval, ci = ci, k = k, m = m, fun = fun, adjust = adjust, statistic = statistic)

   class(res) <- "combp"
   return(res)

}
