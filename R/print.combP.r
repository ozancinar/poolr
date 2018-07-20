print.combP <- function(x, ...) {
  
   cat(" Combined p-value =", round(x$p, 3), "\n",
       "Test Statistic =", round(x$testStat, 3), "\n",
       "Adjustment =", x$adjust, "\n")
  
}
