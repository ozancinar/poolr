print.combp <- function(x, digits=3, ...) {

   if (x$fun %in% c("fisher", "invchisq"))
      info <- paste0("test statistic = ", round(x$statistic, 3), " ~ chi-square(", round(attr(x$statistic, "df"), digits), ")")

   if (x$fun == "stouffer")
      info <- paste0("test statistic = ", round(x$statistic, 3), " ~ N(0,1)")

   if (x$fun %in% c("bonferroni", "tippett"))
      info <- paste0("minimum p-value = ", round(x$statistic, 3))

   if (x$fun == "binotest") {
      if (x$adjust %in% c("nyholt", "liji", "gao", "galwey", "user")) {
         info <- paste0("number of significant tests = ", round(x$statistic * x$m / x$k), " (adjusted based on m)")
      } else {
         info <- paste0("number of significant tests = ", x$statistic)
      }
   }

   if (x$adjust %in% c("nyholt", "liji", "gao", "galwey", "user")) {

      if (x$adjust == "nyholt")
         x$adjust <- "Nyholt, 2004"

      if (x$adjust == "liji")
         x$adjust <- "Li and Ji, 2005"

      if (x$adjust == "gao")
         x$adjust <- "Gao, 2008"

      if (x$adjust == "galwey")
         x$adjust <- "Galwey, 2009"

      if (x$adjust == "user")
         x$adjust <- "user-defined"   

      x$adjust <- paste0("effective number of tests (m) = ", x$m, " (", x$adjust, ")")

   }

   if (x$adjust == "generalized") {

      if (x$fun == "fisher")
         x$adjust <- "Brown's method"
      if (x$fun == "invchisq")
         x$adjust <- "generalized inverse chi-squared method"
      if (x$fun == "stouffer")
         x$adjust <- "Strube's method"

   }

   if (is.null(x$ci)) {
      cat("combined p-value:", round(x$p, digits), "\n")
   } else {
      cat("combined p-value:", round(x$p, digits), paste0("(95% CI: ", round(x$ci[1], digits), ", ", round(x$ci[2], digits), ")"), "\n")
   }

   cat(info, "\n")
   cat("adjustment:", x$adjust, "\n")

   invisible()

}
