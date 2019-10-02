
binotest <- function(p, adjust = "none", m, R, alpha = 0.05, size = 10000, seed, 
                     type = 2, emp.loop = FALSE, emp.step, ...) {

   # checking whether the number of p-values match with the dimensions of R.
   if (!missing(R) && length(p) != nrow(R))
      stop("the number of p-values to be combined does not match with the dimensions of the correlation matrix provided as R.")
   
   # checking if dependence is accounted for when a correlation matrix is supplied.
   if (adjust == "none" && !missing(R))
      warning("although you specified a correlation matrix with argument R, you chose to assume independence with adjust == 'none'. If you want to account for dependence, please specify an adjustment method. See ?binotest for details.")
   
   # if the adjust is not given, it will be set to "none".
   if (missing(adjust)) 
      adjust <- "none"
   
   adjust <- match.arg(adjust, c("none", "nyholt", "liji", "gao", "galwey", "empirical"))
   
   k <- length(p)
   r <- sum(p <= alpha)
   
   # if m is provided by the user, then we don't need to check the adjustment method.
   if (!missing(m)) {
      m <- m
      info <- paste("number of significant individual tests:", round(r * m / k), " (adjusted with meff)")
      # test_stat <- dbinom(round(r * m / k), m, alpha)
      pooled_p <- sum(dbinom(round(r * m / k):m, m, alpha))
      adjust <- paste0("meff = ", m, " (user defined)")
      
      # warning the user if the user-defined m is larger than the number of p-values.
      if (m > k)
         warning("the user-defined effective number of test is larger than the number of p-values that were combined.")
   
   } else {
      if (adjust == "none") {
         info <- paste("number of significant individual tests:", r)
         # test_stat <- dbinom(r, k, alpha)
         pooled_p <- sum(dbinom(r:k, k, alpha))
         adjust <- "none"
      } else if (adjust %in% c("nyholt", "liji", "gao", "galwey")) {
         m <- meff(R = R, method = adjust)
         info <- paste("number of significant individual tests:", round(r * m / k), " (adjusted with meff)")
         # test_stat <- dbinom(round(r * m / k), m, alpha)
         pooled_p <- sum(dbinom(round(r * m / k):m, m, alpha))
         adjust <- paste0("meff = ", m, " (", adjust, ")")
      } else if (adjust == "empirical") {
         # checking if the user wants to use a stepwise algorithm for empirical 
         # adjustment.
         if (!missing(emp.step)) {
            # checking if emp.step is a list.
            if (!is.list(emp.step))
               stop("emp.dist should be a list.")
            
            # checking if there are two vectors in the emp.step.
            if (!length(emp.step) == 2)
               stop("emp.dist should include two lists. Please see ?binotest.")
            
            # checking if the lengths of the vectors in the emp.step are correct.
            if (!(length(emp.step[[1]]) - length(emp.step[[2]])) %in% c(-1, 1))
               stop("the lengths of the vectors in emp.step are not correct. Please see ?binotest.")
            
            # checking if the vectors have names.
            if (is.null(names(emp.step))) {
               # if they don't have names, they will be named.
               names(emp.step)[which.max(unlist(lapply(emp.step, length)))] <- "size"
               names(emp.step)[which.min(unlist(lapply(emp.step, length)))] <- "thres"
            } else {
               # if the names of the vectors are not correct, they will be corrected.
               if (!all(names(emp.step) %in% c("size", "thres"))) {
                  names(emp.step)[which.max(unlist(lapply(emp.step, length)))] <- "size"
                  names(emp.step)[which.min(unlist(lapply(emp.step, length)))] <- "thres"
               }
            }
            # setting the seed before starting to generate sets to obtain mutually exclusive sets.
            if (!missing(seed))
               set.seed(seed)
            
            for (i in 1:length(emp.step$thres)) {
               size <- emp.step$size[i]
               emp_dist <- empirical(R = R, method = "binotest", type = type, 
                                     size = size, emp.loop = emp.loop)
               
               comb_p_tmp <- sum(dbinom(r:k, k, alpha))
               pooled_p_tmp <- (sum(emp_dist <= comb_p_tmp) + 1) / (size + 1)
               
               if (pooled_p_tmp > emp.step$thres[i]) {
                  info <- paste("number of significant individual tests:", r)
                  # test_stat <- comb_p_tmp
                  pooled_p <- pooled_p_tmp
                  ci <- as.numeric(binom.test(sum(emp_dist <= comb_p_tmp) + 1, size + 1)$conf.int)
                  adjust <- "empirical"
                  break
               }
            }
            
            # if the threshold could not be achieved in the loop, then we are going
            # to use the largest sample size. 
            if (!exists("pooled_p")) {
               size <- tail(emp.step$size, n = 1)
               emp_dist <- empirical(R = R, method = "binotest", type = type, 
                                     size = size, seed = seed, emp.loop = emp.loop)
               
               comb_p_tmp <- sum(dbinom(r:k, k, alpha))
               pooled_p <- (sum(emp_dist <= comb_p_tmp) + 1) / (size + 1)
               ci <- as.numeric(binom.test(sum(emp_dist <= pooled_p) + 1, size + 1)$conf.int)
               info <- paste("number of significant individual tests:", r)
               # test_stat <- comb_p_tmp
               adjust <- "empirical"
            }
         } else {
            tmp <- list(...)
            
            # if an empirical distribution is not provided by the user, we will use 
            # empirical() to generate an empirical distribution.
            if (is.null(tmp$emp.dis)) {
               emp_dist <- empirical(R = R, method = "binotest", type = type, 
                                     size = size, seed = seed, emp.loop = emp.loop)
               
            } else { # otherwise, the function will use the user-given empirical distribution.
               emp_dist <- tmp$emp.dist
            }
            
            comb_p_tmp <- sum(dbinom(r:k, k, alpha))
            pooled_p <- (sum(emp_dist <= comb_p_tmp) + 1) / (size + 1)
            ci <- as.numeric(binom.test(sum(emp_dist <= comb_p_tmp) + 1, size + 1)$conf.int)
            info <- paste("number of significant individual tests:", r)
            # test_stat <- comb_p_tmp
            adjust <- "empirical"
         }
      }
   }
   
   if (exists("ci")) {
      res <- list(p = pooled_p, ci = ci, info = info, adjust = adjust)
   } else {
      res <- list(p = pooled_p, info = info, adjust = adjust)
   }
   
   class(res) <- "combp"
   return(res)
   
}
