empirical <- function(p, R, method, type, size = 10000, seed, ...) {

   if (!type %in% c(-1, 1, 2))
      stop("the type of the tests entered is not valid.")

   k <- length(p)

   if (missing(R))
      R <- diag(1, k)

   R <- nearPD(R)$mat

   if (!missing(seed))
      set.seed(seed)

   z <- mvrnorm(size, mu = rep(0, k), Sigma = R)

   if (type == -1) {
      pVals <- pnorm(z, lower.tail = TRUE)
   } else if (type == 1) {
      pVals <- pnorm(z, lower.tail = FALSE)
   } else if (type == 2) {
      pVals <- 2 * pnorm(abs(z), lower.tail = FALSE)
   }

   if (method == "bonferroni")
      emp <- apply(pVals, 1, function(x) bonferroni(x)$testStat)

   if (method == "tippett")
      emp <- apply(pVals, 1, function(x) tippett(x)$testStat)

   if (method == "binotest")
      emp <- apply(pVals, 1, function(x) binotest(x)$testStat)

   if (method == "fisher")
      emp <- apply(pVals, 1, function(x) fisher(x)$testStat)

   if (method == "stouffer")
      emp <- apply(pVals, 1, function(x) stouffer(x)$testStat)

   if (method == "invchisq")
      emp <- apply(pVals, 1, function(x) invchisq(x)$testStat)

   return(emp)

}
