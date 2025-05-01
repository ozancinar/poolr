empirical <- function(R, method, side = 2, size = 10000, batchsize, ...) {

   # check if 'R' is specified
   if (missing(R))
      stop("Argument 'R' must be specified.", call.=FALSE)

   # match 'method' argument
   method <- match.arg(method, c("fisher", "stouffer", "invchisq", "binomtest", "bonferroni", "tippett"))

   # checks for 'side' argument
   .check.side(side)

   # checks for 'R' argument
   R <- .check.R(R, checksym = TRUE, checkna = TRUE, checkpd = TRUE, nearpd = TRUE, checkcor = TRUE, checkdiag = TRUE, isbase = FALSE)

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

   if (batchsize < 1 || batchsize > size)
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

   if (isTRUE(ddd$verbose))
      pbar <- txtProgressBar(min=0, max=length(batchsizes), style=3)

   for (i in seq_along(batchsizes)) {

      if (isTRUE(ddd$verbose))
         setTxtProgressBar(pbar, i)

      if (is.null(ddd$mvnmethod)) {
         mvnmethod <- "mvt_eigen"
      } else {
         mvnmethod <- ddd$mvnmethod
      }

      z <- try(.simmvn(batchsizes[i], Sigma = R, mvnmethod = mvnmethod))

      if (inherits(z, "try-error"))
         stop("Matrix to be generated is too large. Try setting 'batchsize' (or to a lower number if it was set).", call.=FALSE)

      if (side == 1)
         p <- pnorm(z, lower.tail = FALSE)

      if (side == 2)
         p <- 2 * pnorm(abs(z), lower.tail = FALSE)

      emp.dist[(batchpos[i]+1):batchpos[i+1]] <- eval(fcall)

   }

   if (isTRUE(ddd$verbose))
      close(pbar)

   return(emp.dist)

}
