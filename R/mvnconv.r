mvnconv <- function(R, side = 2, target, cov2cor = FALSE) {

   # check if 'R' is specified
   if (missing(R))
      stop("Argument 'R' must be specified.", call.=FALSE)

   # checks for 'R' argument
   if (is.vector(R)) {
      R <- .check.R(R, checksym=FALSE, checkna=FALSE, checkpd=FALSE, checkcor=TRUE, isbase=FALSE)
   } else {
      R <- .check.R(R, checksym=TRUE, checkna=FALSE, checkpd=FALSE, checkcor=TRUE, isbase=FALSE)
   }

   # get name of calling function (character(0) if called from global environment)
   call.fun <- as.character(sys.call(-1)[1])
   call.fun <- gsub("^poolr::", "", call.fun)

   if (isTRUE(call.fun %in% c("fisher", "stouffer", "invchisq"))) {

      if (missing(target)) {

         # for fisher(), stouffer(), and invchisq(), set the default 'target' if it is not specified

         if (call.fun == "fisher")
            target <- "m2lp"
         if (call.fun == "stouffer")
            target <- "z"
         if (call.fun == "invchisq")
            target <- "chisq1"

      }

   } else {

      # when calling mvnconv() from the global environment or some other function besides fisher(), stouffer(), or invchisq(), must specify 'target'

      if (missing(target))
         stop("Argument 'target' must be specified.")

   }

   target <- match.arg(target, c("m2lp", "z", "chisq1", "p"))

   # check for incompatibility between poolr base function and the specified target (only when adjust = "generalized")

   if (length(call.fun) > 0L && call.fun %in% c("fisher", "stouffer", "invchisq")) {
      # need this in case the 'adjust' argument is abbreviated
      call.fun.args <- as.list(match.call(definition = sys.function(-1), call = sys.call(-1), expand.dots = FALSE))
      adjust <- match.arg(call.fun.args$adjust, c("none", "nyholt", "liji", "gao", "galwey", "empirical", "generalized"))
      if (adjust == "generalized" && ((call.fun == "fisher" && target != "m2lp") || (call.fun == "stouffer" && target != "z") || (call.fun == "invchisq" && target != "chisq1")))
         warning(paste0("Using mvnconv(..., target=\"", target, "\") is not compatible with ", call.fun, "()."))
   }

   # checks for 'side' argument

   .check.side(side)

   # set correct column of 'mvnlookup' for converting values in R to target values

   column <- pmatch(target, c("m2lp", "z", "chisq1", "p"))
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

      if (cov2cor) {
         var <- mvnlookup[1,column]
         covs <- covs / var
      }
         

   } else {

      covs <- mvnlookup[match(R, mvnlookup[,1]), column]

      if (cov2cor) {
         var <- mvnlookup[1,column]
         covs <- covs / var
      }

   }

   return(covs)

}
