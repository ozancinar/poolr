
mvnconv <- function(R, side = 2, target, covtocor = FALSE) {

   # check that R is symmetric

   if (!isSymmetric(R))
      stop("R is not symmetric.")

   # checks for 'side' argument

   if (length(side) != 1)
      stop("Argument 'side' must be of length 1.")

   if (!(side %in% c(1,2)))
      stop("Argument 'side' must be either 1 or 2.")

   # get name of calling function (NULL if called from global environment)

   call_fun <- sys.call(-1)[1]

   if (is.null(call_fun)) {

      # when calling from global environment, must specify 'target'
      if (missing(target))
         stop("Argument 'target' must be specified.")

   } else {

      if (missing(target)) {

         call_fun <- as.character(call_fun)

         if (!(call_fun %in% c("fisher", "stouffer", "invchisq")))
            stop("Argument 'target' must be specified.")

         if (call_fun == "fisher")
            target <- "m2lp"
         if (call_fun == "stouffer")
            target <- "z"
         if (call_fun == "invchisq")
            target <- "chisq1"

      }

   }

   column <- pmatch(target, c("m2lp", "z", "chisq1", "p"))

   if (is.na(column))
      stop("Unknown 'target' specified.")

   column <- column * 2

   if (side == 2)
      column <- column + 1

   # get lower triangular part of R
   r <- R[lower.tri(R, diag=TRUE)]

   # round correlations to two decimals
   r <- round(r, 2)

   # replace -1 correlations with -0.99
   r[r == -1] <- -0.99

   # convert correlations into covariances for the chosen target
   covs <- matrix(NA, nrow = nrow(R), ncol = ncol(R))
   covs[lower.tri(covs, diag=TRUE)] <- poolr::mvnlookup[match(r, poolr::mvnlookup[, 1]), column]
   covs[upper.tri(covs)] <- t(covs)[upper.tri(covs)]

   if (covtocor)
      covs <- cov2cor(covs)

   return(covs)

}
