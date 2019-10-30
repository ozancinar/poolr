\name{mvnconv}
\alias{mvnconv}
\title{Converting MVN Correlations to Covariances with Another Distribution}
\description{Function to convert the correlations between the multivariate test statistics to covariances based on a distribution that is assumed by the combining p-values method.}
\usage{
mvnconv(R, side = 2, target, cov2cor = FALSE)
}
\arguments{
   \item{R}{a matrix that includes the correlations among the test statistics that follow a multivariate normal distribution.}
   \item{side}{the type of the test statistics (i.e., the test statistics are obtained from one- or two-sided tests). As default it is set to 2 for two-sided tests. Alternatively, it can bet to 1 for one-sided tests.}
   \item{target}{the target distribution for which the covariances are calculated. It can be set to one of the three possibilities which are "m2lp", "z", or "chisq1" for Brown's method, Stouffer's method, or Chi-Square distribution with 1 df, respectively.}
   \item{cov2cor}{a logical argument to convert the output covariance matrix to correlation matrix.}
}
\details{
   The function can be used to convert the correlations among the test statistics that follow a multivariate normal distribution to covariances assuming a distribution based on the method that is going to be used for combining p-values.
}
\value{
   The function returns the covariance (or, optionally correlation) matrix.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\examples{
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
mvnconv(R, target = "m2lp")
mvnconv(R, target = "m2lp", cov2cor = TRUE)
}
\keyword{htest}
