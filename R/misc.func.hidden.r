.check.p <- function(p) {

   if (missing(p))
      stop("Argument 'p' must be specified.", call.=FALSE)

   if (!is.vector(p) || !is.numeric(p))
      stop("Argument 'p' must be a numeric vector.", call.=FALSE)

   if (any(is.na(p)))
      stop("Values in 'p' vector must not contain NAs.", call.=FALSE)

   if (any(p < 0) || any(p > 1))
      stop("Values in 'p' vector (i.e., the p-values) must be between 0 and 1.", call.=FALSE)

}

.check.R <- function(R, k, adjust, fun) {

   # check that 'R' is a symmetric matrix
   if (!is.matrix(R) || !isSymmetric(unname(R)))
      stop("Argument 'R' must be a symmetric matrix.", call.=FALSE)

   # check if 'R' contains NAs
   if (any(is.na(R)))
      stop("Values in 'R' vector must not contain NAs.", call.=FALSE)

   # check that dimensions of 'R' match the length of 'p'
   if (k != nrow(R))
      stop("Length of 'p' vector (", k, ") does not match the dimensions of the 'R' matrix (", nrow(R), ",", ncol(R), ").", call.=FALSE)

   # check if user specified 'R' argument but no adjustment method
   if (adjust == "none")
      warning("Although argument 'R' was specified, no adjustment method was chosen via the 'adjust' argument.\n  To account for dependence, specify an adjustment method. See help(", fun, ") for details.", call.=FALSE)

}

.check.emp.step <- function(emp.step, size, tmp) {

   # get name of calling function
   call.fun <- as.character(sys.call(-1)[1])

   # if 'emp.step' is not used, set it to the chosen size with threshold 0
   if (missing(emp.step) || !is.null(tmp$emp.dist))
      emp.step <- list(size = size, thres = 0)

   # check if 'emp.step' is a list
   if (!is.list(emp.step))
      stop("Argument 'emp.step' must be a list. See help(", call.fun, ").", call.=FALSE)

   # check if there are two vectors in 'emp.step'
   if (!length(emp.step) == 2)
      stop("Argument 'emp.step' must contain two vectors. See help(", call.fun, ").", call.=FALSE)

   # get lengths of the two vectors in 'emp.step'
   len <- sapply(emp.step, length)

   # check that vectors in 'emp.step' are not of length 0
   if (any(len == 0))
      stop("Vectors in 'emp.step' must have positive lengths. See help(", call.fun, ").", call.=FALSE)

   if (len[1] == len[2]) {

      # if vectors in 'emp.step' are of the same length ...

      if (is.null(names(emp.step))) {

         # if vectors don't have names, check which vector is the threshold vector (all its values must be <= 1)

         is.thres1 <- all(emp.step[1] <= 1)
         is.thres2 <- all(emp.step[2] <= 1)

         if ((is.thres1 && is.thres2) || (!is.thres1 && !is.thres2))
            stop("Cannot identify 'size' and 'thres' vectors in 'emp.step'. See help(", call.fun, ").", call.=FALSE)

         if (is.thres1) {
            names(emp.step)[1] <- "thres"
            names(emp.step)[2] <- "size"
         } else {
            names(emp.step)[1] <- "size"
            names(emp.step)[2] <- "thres"
         }

      } else {

         # if vectors do have names, check that they are 'size' and 'thres'

         if (any(!(names(emp.step) %in% c("size", "thres"))))
            stop("Vectors in 'emp.step' must be named 'size' and 'thres'. See help(", call.fun, ").", call.=FALSE)

      }

      # set last threshold value to 0
      emp.step$thres[length(emp.step$thres)] <- 0

   } else {

      # if vectors in 'emp.step' are NOT of the same length ...

      # check that the lengths of the vectors in 'emp.step' only differ by 1 element

      if (abs(len[1] - len[2]) != 1)
         stop("The lengths of the vectors in 'emp.step' must only differ by one element. See help(", call.fun, ").", call.=FALSE)

      # if so, the longer vector is for 'size', the other for 'thres' (regardless of how they are named)

      if (len[1] > len[2]) {
         names(emp.step)[1] <- "size"
         names(emp.step)[2] <- "thres"
      } else {
         names(emp.step)[1] <- "thres"
         names(emp.step)[2] <- "size"
      }

      # set missing treshold value to 0
      emp.step$thres <- c(emp.step$thres, 0)

   }

   return(emp.step)

}
