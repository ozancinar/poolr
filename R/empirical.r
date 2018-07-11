empirical <- function(R, method, type, size = 10000, seed, ...) {

   # match method argument
   method <- match.arg(method, c("fisher", "stouffer", "invchisq", "binotest", "bonferroni", "tippett"))

   # check type argument
   if (!type %in% c(1, 2))
      stop("Argument 'type' must be either 1 or 2.")

   k <- nrow(R)

   # check that R is symmetric
   if (!isSymmetric(R))
      stop("R is not symmetric.")

   # ensure that the correlation matrix is positive semi-definite
   R <- nearPD(R)$mat

   if (!missing(seed))
      set.seed(seed)

   z <- mvrnorm(size, mu = rep(0, k), Sigma = R)

   if (type == 1) {
      pVals <- pnorm(z, lower.tail = FALSE)
   } else if (type == 2) {
      pVals <- 2 * pnorm(abs(z), lower.tail = FALSE)
   }

   if (method == "fisher")
      emp <- apply(pVals, 1, function(x) fisher(x)$testStat)

   if (method == "stouffer")
      emp <- apply(pVals, 1, function(x) stouffer(x)$testStat)

   if (method == "invchisq")
      emp <- apply(pVals, 1, function(x) invchisq(x)$testStat)

   if (method == "binotest")
      emp <- apply(pVals, 1, function(x) binotest(x)$testStat)

   if (method == "bonferroni")
      emp <- apply(pVals, 1, function(x) bonferroni(x)$testStat)

   if (method == "tippett")
      emp <- apply(pVals, 1, function(x) tippett(x)$testStat)

   return(emp)

}
