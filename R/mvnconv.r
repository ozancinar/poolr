mvnconv <- function(R, side = 2, target, cov2cor = FALSE) {

   # check if 'R' is specified
   if (missing(R))
      stop("Argument 'R' must be specified.", call.=FALSE)

   # get name of calling function (NULL if called from global environment)
   call.fun <- sys.call(-1)[1]

   if (is.null(call.fun)) {

      # when calling mvnconv() from global environment, must specify 'target'

      if (missing(target))
         stop("Argument 'target' must be specified.")

   } else {

      # if mvnconv() is called inside one of the poolr functions, check that 'R' is a symmetric matrix

      call.fun <- as.character(call.fun)

      if (call.fun %in% c("fisher", "invchisq", "stouffer", "bonferroni", "tippett", "binotest")) {
         if (!is.matrix(R) || !isSymmetric(unname(R)))
            stop("Argument 'R' must be a symmetric matrix.")
      }

      if (missing(target)) {

         if (!(call.fun %in% c("fisher", "stouffer", "invchisq")))
            stop("Argument 'target' must be specified.")

         if (call.fun == "fisher")
            target <- "m2lp"
         if (call.fun == "stouffer")
            target <- "z"
         if (call.fun == "invchisq")
            target <- "chisq1"

      }

   }

   # checks for 'side' argument

   if (length(side) != 1)
      stop("Argument 'side' must be of length 1.")

   if (!(side %in% c(1,2)))
      stop("Argument 'side' must be either 1 or 2.")

   column <- pmatch(target, c("m2lp", "z", "chisq1", "p"))

   if (is.na(column))
      stop("Unknown 'target' specified.")

   column <- column * 2

   if (side == 2)
      column <- column + 1

   # round elements in 'R' to two decimals (since mvnlookup[,1] values are in .01 steps)
   R <- round(R, 2)

   # replace -1 elements in 'R' with -0.99
   R[R == -1] <- -0.99

   mvnlookup <- get(data(mvnlookup, package="poolr", envir = environment()))

   if (is.matrix(R)) {

      # get lower triangular part of R
      r <- R[lower.tri(R, diag=TRUE)]

      # convert correlations to covariances for the chosen target
      covs <- matrix(NA, nrow = nrow(R), ncol = ncol(R))
      covs[lower.tri(covs, diag=TRUE)] <- mvnlookup[match(r, mvnlookup[,1]), column]
      covs[upper.tri(covs)] <- t(covs)[upper.tri(covs)]

      if (cov2cor)
         covs <- stats::cov2cor(covs)

   } else {

      covs <- mvnlookup[match(R, mvnlookup[,1]), column]

      if (cov2cor) {
         var <- tail(mvnlookup[,column], 1)
         covs <- covs / var
      }

   }

   return(covs)

}
