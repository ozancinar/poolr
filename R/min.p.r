
min.p <- function(p, adjust = "none") {
   pooled.p <- min(p)
   res <- list(p = pooled.p, adjust = paste0(adjust, " "))
   return(res)
}