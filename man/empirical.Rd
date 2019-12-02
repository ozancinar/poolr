\name{empirical}
\alias{empirical}
\title{Generate Empirically-Derived Null Distributions}
\description{Function to derive the empirical distribution of the combined p-value under the null hypothesis.}
\usage{
empirical(R, method, side, size = 10000, emploop = FALSE, \dots)
}
\arguments{
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values. If it is NULL then the p-values are assumed to be independent.}
   \item{method}{an argument to specify the method to be used to derive the empirical distribution of the p-values under the null hypothesis. It can be set to one of the base function name, "fisher", "stouffer", "invchisq", "binotest", "bonferroni", or "tippett".}
   \item{side}{the side of the hypothesis test. It can take values 1 or 2 for one-sided or two-sided p-values, respectively.}
   \item{size}{size of the empirical distribution. It is set to 10000 as default.}
   \item{emploop}{a logical argument to decide to generate the empirical distribution in a loop. If the size of the empirical distribution to be generated is too large, we recommend you to use emploop = TRUE.}
   \item{\dots}{other arguments.}
}
\details{
   This function can be used to derive the empirical distribution of the combined p-value under the null hypothesis. The function can be used with one of the six base functions in poolR by stating the name of the function with the "method" argument.
}
\value{
   The function returns a vector that includes a set of combined p-values that represents the empirical distribution under the null hypothesis.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Mehta, C. R., & Patel, N. R. (1998). \emph{Exact inference for categorical data}. Encyclopedia of biostatistics, 2, 1411-1422.
}
\examples{
# create an example correlation matrix with constant (positive) correlations
R <- matrix(0.3, nrow = 4, ncol = 4)
diag(R) <- 1

# generate empirical distributions for fisher() with different features
emp_tmp <- empirical(R, method = "fisher", side = 1)  # for one-sided p-values
emp_tmp <- empirical(R, method = "fisher", side = 2)  # for two-sided p-values
emp_tmp <- empirical(R, method = "fisher", side = 2, emploop = TRUE)  # using a loop
}
\keyword{htest}
