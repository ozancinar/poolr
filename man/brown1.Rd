\name{brown1}
\alias{brown1}
\title{Brown's Method}
\description{Function to get Brown's correlation matrix for a given correlation matrix for one-sided p-values.}
\usage{
brown1(R, \dots)
}
\arguments{
   \item{R}{a symmetric matrix that includes pairwise correlations between one-sided p-values.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to apply Brown's method to adjust Fisher's method. The function itself can be used to calculate a combined p-value under the dependence of the individual one-sided p-values. Moreover, this function is also embedded inside the fisher() to adjust that function to adjust the Fisher's method with Brown's method.
}
\value{
   The function returns a matrix that includes correlation values with respect to Brown's method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Brown, MB (1975): \emph{400: A method for combining non-independent, one-sided tests of significance}. Biometrics, 987-992.
}
\examples{
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
brown1(R)
}
\keyword{htest}
