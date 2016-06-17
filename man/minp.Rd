\name{minp}
\alias{minp}
\title{Minimum p}
\description{Function to combine the p-values with the strongest evidence, i.e. the smallest p-values.}
\usage{
minp(p, adjust, pca.method, R, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. As default, the Bonferroni correction is applied to the minimum p-value. Optionally, number of significant tests can be used instead of the Bonferroni correction.}
   \item{pca.method}{the method to be used to calculate the number of effective tests.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to combine a set of p-values with the strongest evidence which is determined by the smallest p-value.
}
\value{
   The function returns the (pooled) p-value based on minimum p-value method.
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
minp(p = p)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
minp(p = p, adjust = "m.eff", pca.method = "nyholt", R = R)
}
\keyword{htest}
