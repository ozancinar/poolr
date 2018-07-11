binotest <- function(p, adjust = "none", m, R, alpha = 0.05, size = 10000, seed, type = 2, ...) {

   k <- length(p)
   r <- sum(p <= alpha)

   testStat <- dbinom(r, k, alpha)

   if (adjust == "none") {

      pooled.p <- sum(dbinom(r:k, k, alpha))
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

      testStat <- dbinom(round(r * m / k), m, alpha)
      pooled.p <- sum(dbinom(round(r * m / k):m, m, alpha))

   }

   if (adjust == "empirical") {

      tmp <- list(...)
      if (is.null(tmp$emp.dis)) {
         emp.dist <- empirical(R = R, method = "binotest", type = type, size = size, seed = seed)
      } else {
         emp.dist <- tmp$emp.dist
      }

      pooled.p <- sum(emp.dist <= testStat) / length(emp.dist)
      adjust <- "empirical"

   }

   res <- list(p = pooled.p, testStat = testStat, adjust = adjust)
   class(res) <- "combP"
   return(res)

}
