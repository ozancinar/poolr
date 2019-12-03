empirical <- function(R, method, side = 2, size = 10000, batchsize, ...) {

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
   if (any(eigen(R)$values <= 0)) {
      R <- as.matrix(Matrix::nearPD(R, corr=TRUE)$mat)
      warning("Matrix 'R' is not positive definite. Used Matrix::nearPD() to make 'R' positive definite.", call.=FALSE)
   }

   ddd <- list(...)

   if (is.null(ddd$alpha)) {
      alpha <- .05
   } else {
      alpha <- ddd$alpha
   }

   k <- nrow(R)
   mu <- rep(0, k)

   if (missing(batchsize) || is.null(batchsize))
      batchsize <- size

   if (batchsize < 1 && batchsize > size)
      stop("Argument 'batchsize' must be between 1 and the value of the 'size' argument.")

   emp.dist <- rep(NA_real_, size)

   fcall <- parse(text=paste0("apply(p, 1, function(x) .", method, "(x, k, alpha))"))

   if (size %% batchsize == 0) {
      batches <- size / batchsize
      batchsizes <- rep(batchsize, batches)
   } else {
      batches <- floor(size / batchsize) + 1
      batchsizes <- c(rep(batchsize, batches - 1), size %% batchsize)
   }

   batchpos <- c(0, cumsum(batchsizes))

   #return(list(batches=batches, batchsize=batchsize, batchsizes=batchsizes, batchpos=batchpos))

   if (!is.null(ddd$verbose) && ddd$verbose)
      pbar <- txtProgressBar(min=0, max=length(batchsizes), style=3)

   for (i in seq_along(batchsizes)) {

      if (!is.null(ddd$verbose) && ddd$verbose)
         setTxtProgressBar(pbar, i)

      z <- try(.simmvn(batchsizes[i], mu = mu, Sigma = R))

      if (inherits(z, "try-error"))
         stop("Matrix to be generated is too large. Try setting 'batchsize' (or to a lower number if it was set).", call.=FALSE)

      if (side == 1)
         p <- pnorm(z, lower.tail = FALSE)

      if (side == 2)
         p <- 2 * pnorm(abs(z), lower.tail = FALSE)

      emp.dist[(batchpos[i]+1):batchpos[i+1]] <- eval(fcall)

   }

   if (!is.null(ddd$verbose) && ddd$verbose)
      close(pbar)

   return(emp.dist)

}
