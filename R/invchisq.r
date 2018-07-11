invchisq <- function(p, adjust = "none", m, R, size = 10000, seed, type = 2, ...) {

   k <- length(p)

   testStat <- sum(qchisq(p, df = 1, lower.tail = FALSE))

   if (adjust == "none") {

      pooled.p <- pchisq(testStat, df = k, lower.tail = FALSE)
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
      pooled.p <- pchisq(testStat, df = m, lower.tail = FALSE)

   }

   if (adjust == "empirical") {

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = "invchisq", type = type, size = size, seed = seed)
      } else {
         emp.dist <- tmp$emp.dist
      }

      pooled.p <- sum(emp.dist >= testStat) / length(emp.dist)
      adjust <- "empirical"

   }

   res <- list(p = pooled.p, testStat = testStat, adjust = adjust)
   class(res) <- "combP"
   return(res)

}
