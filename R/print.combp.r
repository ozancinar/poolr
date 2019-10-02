
print.combp <- function (x, ...) {

   if (length(x) == 3) {
      cat(" combined p-value:", round(x$p, 3), "\n", 
          "adjustment:", x$adjust, "\n",
          "info:", x$info,  
          "\n")
   } else if (length(x) == 4) {
      cat(" combined p-value:", round(x$p, 3), paste0(" (", round(x$ci[1], 3), ", ", round(x$ci[2], 3), ")"), "\n", 
          "adjustment:", x$adjust, "\n",
          "info:", x$info, 
          "\n")
   }
   
}
