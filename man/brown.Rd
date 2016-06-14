\name{brown}
\alias{brown}
\title{Brown's Method}
\description{Function to carry out Brown's method to adjust Fisher's method to consider the dependencies between the p-values.}
\usage{
brown(p, R, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to apply Brown's method to adjust Fisher's method. The function itself can be used to calculate a combined p-value under the dependence of the individual p-values. Moreover, this function is also embedded inside the fisher() to adjust that function to adjust the Fisher's method with Brown's method.
}
\value{
   The function returns the (pooled) p-value based on Brown's method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Brown, MB (1975): \emph{400: A method for combining non-independent, one-sided tests of significance}. Biometrics, 987-992.
}
\examples{
p <- c(.013, .022, .124, .278)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
brown(p = p, R = R)
}
\keyword{htest}
