empirical <- function(R, method, side, size = 10000, emploop = FALSE, ...) {

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

   mu <- rep(0, nrow(R))

   if (!emploop) {

      z <- MASS::mvrnorm(size, mu = mu, Sigma = R)

      if (side == 1)
         p <- pnorm(z, lower.tail = FALSE)

      if (side == 2)
         p <- 2 * pnorm(abs(z), lower.tail = FALSE)

      #emp.dist <- apply(p, 1, function(x) do.call(method, list(x))$p)
      emp.dist <- eval(parse(text=paste0("apply(p, 1, function(x) ", method, "(x)$p)")))

   } else {

      emp.dist <- rep(NA_real_, size)

      fcall <- parse(text=paste0(method, "(p)$p"))

      for (i in seq_len(size)) {

         z <- MASS::mvrnorm(1, mu = mu, Sigma = R)

         if (side == 1)
            p <- pnorm(z, lower.tail = FALSE)

         if (side == 2)
            p <- 2 * pnorm(abs(z), lower.tail = FALSE)

         emp.dist[i] <- eval(fcall)

      }

   }

   return(emp.dist)

}
