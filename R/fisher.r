fisher <- function(p, adjust = "none", m, R, size = 10000, seed, type = 2, ...) {

   if(!adjust %in% c("none", "nyholt", "liji", "gao", "galwey", "empirical", "brown", "strube") & missing(m)) {
     stop("adjustment method is not correct and the user-defined effective number of tests is missing.")
   }


   k <- length(p)

   testStat <- -2 * sum(log(p))

   if (adjust == "none") {

      pooled.p <- pchisq(testStat, df = 2*k, lower.tail = FALSE)
      adjust <- "none"

   }

   if (adjust %in% c("nyholt", "liji", "gao", "galwey")) {

      if (!missing(m)) {
         m <- m
         adjust <- paste0(m, " (user defined)")
      } else {
         m <- meff(R = R, method = adjust)
         adjust <- paste0("meff (", adjust, ")")
      }

      testStat <- testStat * (m / k)
      pooled.p <- pchisq(testStat, df = 2 * m, lower.tail = FALSE)

   }

   if (adjust == "brown") {

      tmp <- list(...)
      if (is.null(tmp$brownCov)) {
         covs <- brown(R)
      } else {
         covs <- tmp$brownCov
      }

      expx2 <- 2 * k
      varx2 <- sum(covs)
      fval <- 2 * expx2^2 / varx2
      cval <- varx2 / (2 * expx2)

      testStat <- testStat/cval
      pooled.p <- pchisq(testStat, df = fval, lower.tail = FALSE)
      adjust <- "brown"

   }

   if (adjust == "empirical") {

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = "fisher", type = type, size = size, seed = seed)
      } else {
         emp.dist <- tmp$emp.dist
      }

      pooled.p <- (sum(emp.dist >= testStat) + 1) / (size + 1)
      adjust <- "empirical"

   }

   res <- list(p = pooled.p, testStat = testStat, adjust = adjust)
   class(res) <- "combP"
   return(res)

}
