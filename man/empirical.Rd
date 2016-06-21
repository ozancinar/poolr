\name{empirical}
\alias{empirical}
\title{Emprically-Derived Null Distributions}
\description{Function to derive the empirical distribution of the combined p-value under the null hypothesis.}
\usage{
empirical(p, R, method, size, seed, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values. If it is NULL then the p-values are assumed to be independent.}
   \item{method}{an argument to specify the method to be used to derive the empirical distribution of the p-values under the null hypothesis. This function can be used to adjust one of minimum p-value, binomial test, Fisher's method, or Stouffer's method. Therefore, this argument can be set to "minp", "binotest", "fisher", or "stouffer".}
   \item{size}{size of the empirical distribution. It is set to 10000 as default.}
   \item{seed}{an argument that can be used to set the seed for random number generation.}
   \item{\dots}{other arguments.}
}
\details{
   This function can be used to derive the empirical distribution of the combined p-value under the null hypothesis. The function can be used with minimum p-value, binomial test, Fisher's method, or Stouffer's method. 
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
p <- c(.013, .022, .124, .278)
empirical(p = p, method = "fisher")
p <- c(.013, .022, .124, .278)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
fisher(p = p, R = R, method = "stouffer", size = 1000, seed = 123)
}
\keyword{htest}
