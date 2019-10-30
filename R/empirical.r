empirical <- function(R, method, side, size = 10000, seed, emp.loop = FALSE, ...) {

   # check if 'R' is specified
   if (missing(R))
      stop("Argument 'R' must be specified.", call.=FALSE)

   # match 'method' argument
   method <- match.arg(method, c("fisher", "stouffer", "invchisq", "binotest", "bonferroni", "tippett"))

   # check 'side' argument
   if (!side %in% c(1, 2))
      stop("Argument 'side' must be either 1 or 2.", call.=FALSE)

   # check that 'R' is a symmetric matrix
   if (!is.matrix(R) || !isSymmetric(unname(R)))
      stop("Argument 'R' must be a symmetric matrix.", call.=FALSE)

   # check if 'R' is positive definite; if not, make it
   if (!all(eigen(R)$values > 0)) {
      R <- as.matrix(Matrix::nearPD(R, corr=TRUE)$mat)
      warning("Matrix 'R' is not positive definite. Used Matrix::nearPD() to make 'R' positive definite.", call.=FALSE)
   }

   # set seed if specified
   if (!missing(seed))
      set.seed(seed)

   mu <- rep(0, nrow(R))

   if (!emp.loop) {

      z <- mvrnorm(size, mu = mu, Sigma = R)

      if (side == 1)
         p <- pnorm(z, lower.tail = FALSE)

      if (side == 2)
         p <- 2 * pnorm(abs(z), lower.tail = FALSE)

      if (method == "fisher")
         emp.dist <- apply(p, 1, function(x) fisher(x)$p)

      if (method == "stouffer")
         emp.dist <- apply(p, 1, function(x) stouffer(x)$p)

      if (method == "invchisq")
         emp.dist <- apply(p, 1, function(x) invchisq(x)$p)

      if (method == "binotest")
         emp.dist <- apply(p, 1, function(x) binotest(x)$p)

      if (method == "bonferroni")
         emp.dist <- apply(p, 1, function(x) bonferroni(x)$p)

      if (method == "tippett")
         emp.dist <- apply(p, 1, function(x) tippett(x)$p)

   } else {

      emp.dist <- rep(NA_real_, size)

      for (i in 1:size) {

         z <- mvrnorm(1, mu = mu, Sigma = R)

         if (side == 1)
            p <- pnorm(z, lower.tail = FALSE)

         if (side == 2)
            p <- 2 * pnorm(abs(z), lower.tail = FALSE)

         if (method == "fisher")
            emp.dist[i] <- fisher(p)$p

         if (method == "stouffer")
            emp.dist[i] <- stouffer(p)$p

         if (method == "invchisq")
            emp.dist[i] <- invchisq(p)$p

         if (method == "binotest")
            emp.dist[i] <- binotest(p)$p

         if (method == "bonferroni")
            emp.dist[i] <- bonferroni(p)$p

         if (method == "tippett")
            emp.dist[i] <- tippett(p)$p

      }

   }

   return(emp.dist)

}
