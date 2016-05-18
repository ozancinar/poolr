fisher <- function(p, ...) {
   k <- length(p)
   pchisq(-2*sum(log(p)), df=2*k, lower.tail=FALSE)
}
