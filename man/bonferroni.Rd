\name{bonferroni}
\alias{bonferroni}
\title{Bonferroni Correction}
\description{Function to combine the p-values by using the Bonferroni correction.}
\usage{
bonferroni(p, adjust = "none", pca.method, R,
           size = 10000, seed, type = 2, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. As default, the Bonferroni correction is applied to the minimum p-value. Optionally, number of significant tests can be used instead of the Bonferroni correction.}
   \item{pca.method}{the method to be used to calculate the number of effective tests.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{size}{size of the emprical distribution if the empirically-derived null hypothesis is going to be used.}
   \item{seed}{seed for the random number generation if the empirically-derived null hypothesis is going to be used.}
   \item{type}{type of the p-values used in the tests. As default, it is set to 2 for two-sided tests.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to combine a set of p-values by using the Bonferroni correction and the strongest evidence, i.e. the smallest p-value.
}
\value{
   The function returns the (pooled) p-value based on the Bonferroni correction.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Bland, J. M., & Altman, D. G. (1995). \emph{Multiple significance tests: The Bonferroni method}. British Medical Journal, 310(6973), 170.
}
\examples{
p <- c(.013, .022, .124, .278)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
bonferroni(p)
bonferroni(p, adjust = "m.eff", pca.method = "nyholt", R = R)
}
\keyword{htest}
