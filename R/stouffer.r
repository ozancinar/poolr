stouffer <- function(p, adjust = "none", m, R, size = 10000, seed, type = 2, ...) {

   k <- length(p)

   testStat <- sum(qnorm(p, lower.tail = FALSE)) / sqrt(k)

   if (adjust == "none") {

      pooled.p <- pnorm(testStat, lower.tail = FALSE)
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

      testStat <- testStat * sqrt(m / k)
      pooled.p <- pnorm(testStat, lower.tail = FALSE)

   }

   if (adjust == "strube") {

      testStat <- sum(qnorm(p, lower.tail = FALSE)) / sqrt(sum(R))
      pooled.p <- pnorm(testStat, lower.tail = FALSE)
      adjust <- "strube"

   }

   if (adjust == "empirical") {

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = "stouffer", type = type, size = size, seed = seed)
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
