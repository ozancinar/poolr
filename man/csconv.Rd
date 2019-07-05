\name{csconv}
\alias{csconv}
\title{Converting MVN Correlations to Covariances with Another Distribution Assuming Compound Symmetry}
\description{Function to convert the correlations between the multivariate test statistics to covariances based on a distribution assuming compound symmetry for the correlations among the tests.}
\usage{
csconv(p, target)
}
\arguments{
   \item{p}{a vector of $p$-values that are going to be used to estimate the correaltions among the tests.}
   \item{target}{the target distribution for which the covariances are calculated. It can be set to one of the three possibilities which are "m2lp", "z", or "chisq1" for Brown's method, Stouffer's method, or Chi-Square distribution with 1 df, respectively.}
}
\details{
   The function can be used to estimate the correlations among the tests assuming compound symmetry.
}
\value{
   The function returns the covariance matrix.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\examples{
rho <- 0.5
k <- 10
r <- matrix(rho, k, k); diag(r) <- 1

set.seed(123)
z <- mvrnorm(1, mu = rep(0, k), Sigma = r)
p <- pnorm(z)

csconv(p = p, target = "m2lp")
csconv(p = p, target = "z")
csconv(p = p, target = "chisq1")
}
\keyword{htest}
