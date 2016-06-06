\name{binotest}
\alias{binotest}
\title{Binomial Test}
\description{Function to carry out Binomial Test to combine the p-values.}
\usage{
binotest(p, adjust, pca.method, R, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. Can be "none", as default, to use the default binomial test. Or the test can be adjusted with PCA-based methods by calculating the number of effective tests.}
   \item{pca.method}{the method to be used to calculate the number of effective tests.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to use Binomial dsitribution to combine the p-values. The user can choose to use default parameters for the binomial test depending on the number of p-values to be combined by leaving "adjust" argument as default. Or, the user can choose to adjust the 'n' parameter of the binomial distribution with PCA-based methods to consider the dependencies between the tests.
}
\value{
   The function returns the (pooled) p-value calculated with binomial distribution and the selected adjustment method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Wilkinson, B. (1951). \emph{A statistical consideration in psychological research}. Psychological Bulletin, 48(2), 156-158.
}
\examples{
p <- c(.013, .022, .124, .278)
binotest(p)
p <- c(.013, .022, .124, .278)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
binotest(p, adjust = "m.eff", pca.method = "nyholt", R = R)
}
\keyword{htest}
