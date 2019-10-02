
mvnconv <- function(R, side = 2, target, covtocor = FALSE) {

   # dimension checks
   if (!isSymmetric(R))
      stop("R is not symmetric.")
   
   test_fun <- NULL
   test_fun <- try(sys.call(-1))
   
   if (is.null(test_fun)) {
      if (missing(target)) {
         stop("target is missing")
      } else {
         target <- target
      }
   } else {
      if (missing(target)) {
         fun_info <- as.character(test_fun)
         fun_name <- fun_info[which(fun_info %in% c("fisher", "stouffer", "invchisq"))]
         
         if (fun_name == "fisher") {
            target <- "m2lp"
         } else if (fun_name == "stouffer") {
            target <- "z"
         } else if (fun_name == "invchisq") {
            target <- "chisq1"
         } else {
            stop("mvnconv() was called from a function that is not provided by poolR. This is not available for now.")
         }
      } else {
         target <- target
      }
   }
   
   look <- which(c("m2lp", "z", "chisq1", "p") %in% target)
   look <- look * 2
   if (side == 2) {
      look <- look + 1
   }
   
   data(mvnlookup)
   
   # lower triangular part of R
   r <- R[lower.tri(R, diag=TRUE)]
   
   # round correlations to two decimals
   r <- round(r, 2)
   
   # replace -1 correlations with -.999
   r[r == -1] <- -0.99
   
   # convert correlations into covariances
   covs <- matrix(NA, nrow = nrow(R), ncol = ncol(R))
   covs[lower.tri(covs, diag=TRUE)] <- mvnlookup[match(r, mvnlookup[, 1]), look]
   covs[upper.tri(covs)] <- t(covs)[upper.tri(covs)]
   
   if (covtocor)
      covs <- cov2cor(covs)
   
   return(covs)
   
}
