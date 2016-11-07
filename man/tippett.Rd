\name{tippett}
\alias{tippett}
\title{Tippett's Minimum P-Value Method}
\description{Function to combine the p-values by using Tippett's minimum p-value method.}
\usage{
tippett(p, adjust, pca.method, R, size, seed, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. As default, Tippett's minumum p-value method is applied to the minimum p-value. Optionally, number of effective tests or empirically-derived null distributions can be used instead of the Tippett's minimum p-value method.}
   \item{pca.method}{the method to be used to calculate the number of effective tests.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{size}{size of the emprical distribution if the empirically-derived null hypothesis is going to be used.}
   \item{seed}{seed for the random number generation if the empirically-derived null hypothesis is going to be used.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to combine a set of p-values by using Tippett's minimum p-value method and possible adjustments of it to consider the dependency among the p-values.
}
\value{
   The function returns the (pooled) p-value based on the Tippett's minimum p-value method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Tippett, L. H. C. (1931). \emph{The Methods of Statistics}. The methods of statistics.
}
\examples{
p <- c(.013, .022, .124, .278)
tippett(p = p)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
tippett(p = p, adjust = "m.eff", pca.method = "nyholt", R = R)
}
\keyword{htest}
