print.poolr <- function(x, digits=3, ...) {

   cat("combined p-values with:      ")

   if (x$fun == "fisher")
      cat("Fisher's method\n")
   if (x$fun == "stouffer")
      cat("Stouffer's method\n")
   if (x$fun == "invchisq")
      cat("inverse chi-square method\n")
   if (x$fun == "bonferroni")
      cat("Bonferroni method\n")
   if (x$fun == "tippett")
      cat("Tippett's method\n")
   if (x$fun == "binomtest")
      cat("binomial test\n")

   cat("number of p-values combined:", x$k, "\n")

   if (x$fun %in% c("fisher", "invchisq"))
      testinfo <- paste0("test statistic:              ", round(x$statistic, digits), " ~ chi-square(df = ", round(attr(x$statistic, "df"), digits), ")")

   if (x$fun == "stouffer")
      testinfo <- paste0("test statistic:              ", round(x$statistic, digits), " ~ N(0,1)")

   if (x$fun %in% c("bonferroni", "tippett"))
      testinfo <- paste0("minimum p-value:             ", round(x$statistic, digits))

   if (x$fun == "binomtest") {
      if (x$adjust %in% c("nyholt", "liji", "gao", "galwey", "chen", "user")) {
         testinfo <- paste0("number of significant tests: ", round(x$statistic * x$m / x$k), " (adjusted based on m; at alpha = ", x$alpha, ")")
      } else {
         testinfo <- paste0("number of significant tests: ", x$statistic, " (at alpha = ", x$alpha, ")")
      }
   }

   cat(testinfo, "\n")

   if (x$adjust %in% c("nyholt", "liji", "gao", "galwey", "chen", "user")) {

      if (x$adjust == "nyholt")
         x$adjust <- "Nyholt, 2004"

      if (x$adjust == "liji")
         x$adjust <- "Li & Ji, 2005"

      if (x$adjust == "gao")
         x$adjust <- "Gao, 2008"

      if (x$adjust == "galwey")
         x$adjust <- "Galwey, 2009"

      if (x$adjust == "chen")
         x$adjust <- "Chen & Liu, 2011"

      if (x$adjust == "user")
         x$adjust <- "user-defined"

      x$adjust <- paste0("effective number of tests (m = ", x$m, "; ", x$adjust, ")")
      x$p <- format.pval(x$p, digits)

   }

   if (x$adjust == "generalized") {

      if (x$fun == "fisher")
         x$adjust <- "Brown's method"
      if (x$fun == "invchisq")
         x$adjust <- "Satterthwaite approximation"
      if (x$fun == "stouffer")
         x$adjust <- "Strube's method"

      x$p <- format.pval(x$p, digits)

   }

   if (x$adjust == "empirical") {

      x$adjust <- paste0("empirical distribution (size = ", as.integer(x$size), ")")
      x$p <- paste0(format.pval(x$p, digits), " (95% CI: ", format.pval(x$ci[1], digits), ", ", format.pval(x$ci[2], digits), ")")

   }

   cat("adjustment:                 ", x$adjust, "\n")
   cat("combined p-value:           ", x$p, "\n")

   invisible()

}
