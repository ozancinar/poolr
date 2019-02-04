mvnconv <- function(R, side = 2, target, covtocor = FALSE) {

  # test if the function is called directly or inside another function
  testFun <- try(sys.call(-1))

  if (is.null(testFun)) {
    # if the upper test is true then this function is directly
    if (missing(target)) {
      # if the target is not defined either, then this function cannot work by itself
      stop("Please specify a target distribution. See ?mvnconv for details and examples.")
    } else {
      if (!target %in% c("m2lp", "z", "chisq1", "p")) {
        stop("target should be one of 'm2lp', 'z', 'chisq1', or 'p'.")
      } else {
        look <- which(c("m2lp", "z", "chisq1", "p") %in% target)
      }
    }

  } else {
    # if else, then this function is called inside another function
    if (missing(target)) {
      # if this is the case, this function is called inside a base function and
      # we need to use the name of this base function to define the target

      nameFun <- as.character(testFun)[1] # get name of function from which this
                                          # function is being called

      if (nameFun %in% c("fisher", "stouffer", "invchisq")) {
        look <- which(c("fisher", "stouffer", "invchisq") %in% nameFun)
      } else {
        stop("You are calling this function from a non-poolR function. You can't do this for now.")
      }

    } else {
      if (!target %in% c("m2lp", "z", "chisq1", "p")) {
        stop("target should be one of 'm2lp', 'z', 'chisq1', or 'p'.")
      } else {
        look <- which(c("m2lp", "z", "chisq1", "p") %in% target)
      }
    }
  }

  data(mvnlookup)

  if (side == 1) {
    if (look == 1)
      colNum <- 2
    else if (look == 2)
      colNum <- 4
    else if (look == 3)
      colNum <- 6
    else if (look == 4)
      colNum <- 8
  } else if (side == 2) {
    if (look == 1)
      colNum <- 3
    else if (look == 2)
      colNum <- 5
    else if (look == 3)
      colNum <- 7
    else if (look == 4)
      colNum <- 9
  }

  # dimension checks
  if (!isSymmetric(R))
    stop("R is not symmetric.")

  # lower triangular part of R
  r <- R[lower.tri(R, diag=TRUE)]

  # round correlations to two decimals
  r <- round(r, 3)

  # replace -1 correlations with -.999
  r[r == -1] <- -0.999

  # convert correlations into covariances
  covs <- matrix(NA, nrow = nrow(R), ncol = ncol(R))
  covs[lower.tri(covs, diag=TRUE)] <- mvnlookup[match(r, mvnlookup[, 1]), colNum]
  covs[upper.tri(covs)] <- t(covs)[upper.tri(covs)]

  if (covtocor)
    covs <- cov2cor(covs)

  # return cov matrix
  return(covs)

}
