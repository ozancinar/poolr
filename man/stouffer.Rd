\name{stouffer}
\alias{stouffer}
\title{Stouffer's Method}
\description{Function to carry out Stouffer's method.}
\usage{
Stouffer(p, adjust, pca.method, R, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. Can be "none", as default, to use the default Stouffer's method. Or Stouffer's method can be adjusted with PCA-based methods by calculating the number of effective tests.}
   \item{pca.method}{the method to be used to calculate the number of effective tests.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to apply Stouffer's method. The user can choose to use default Stouffer's method by leaving "adjust" argument as default. Or, the user can choose to adjust the Stouffer's method with PCA-based methods to consider the dependencies between the tests.
}
\value{
   The function returns the (pooled) p-value based on Stouffer's method and the selected adjustment method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Strube, M. J. (1985). \emph{Combining and comparing significance levels from nonindependent hypothesis tests}. Psychological bulletin, 97(2), 334.
}
\examples{
p <- c(.013, .022, .124, .278)
stouffer(p)
p <- c(.013, .022, .124, .278)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
stouffer(p, adjust = "m.eff", pca.method = "nyholt", R = R)
}
\keyword{htest}
