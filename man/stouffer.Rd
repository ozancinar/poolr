\name{stouffer}
\alias{stouffer}
\title{Stouffer's Method}
\description{Function to carry out Stouffer's method.}
\usage{
stouffer(p, adjust = "none", pca.method, R,
         size = 10000, seed, type = 2, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. Can be "none", as default, to use the default Stouffer's method. Or Stouffer's method can be adjusted with PCA-based methods by calculating the number of effective tests. Finally, this argument can be set to "strube" in order to calculate the test statistic under dependence with Strube's method.}
   \item{pca.method}{the method to be used to calculate the number of effective tests.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{size}{size of the emprical distribution if the empirically-derived null hypothesis is going to be used.}
   \item{seed}{seed for the random number generation if the empirically-derived null hypothesis is going to be used.}
   \item{type}{type of the p-values used in the tests. As default, it is set to 2 for two-sided tests.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to apply Stouffer's method. The user can choose to use default Stouffer's method by leaving "adjust" argument as default. Or, the user can choose to adjust the Stouffer's method with PCA-based methods to consider the dependencies between the tests. As a final modification, Strube's method can be used to calculate the test statistic under dependence by setting the "adjust" argument to "strube".
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
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
stouffer(p)
stouffer(p, adjust = "m.eff", pca.method = "nyholt", R = R)
}
\keyword{htest}
