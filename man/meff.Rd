\name{meff}
\alias{meff}
\title{Number of Effective Tests}
\description{Function to calculate the number of effective tests based on one of the four methods.}
\usage{
meff(R, eigen = FALSE, method, \dots)
}
\arguments{
   \item{R}{main input for the function. It can be a matrix for which the number of effective tests is calculated. Or, it can be pre-calculated eigenvalues of the matrix. If this is the case, the second argument eigen should be TRUE.}
   \item{eigen}{states if the main input is the eigen values. If this is true the function will not calculate the eigen values, and use the given values as the eigen values.}
   \item{method}{method to be used to calculate the number of effective tests. Can be one of "nyholt", "liji", "gao", or "galwey".}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to calculate the number of effective tests based on one of the four methods. Du to the problems arise with numerical calculations, in Li & Ji method, a small value is added to the eigen values.
}
\value{
   The function returns the number of effective tests.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Gao, X., Starmer, J., & Martin, E. R. (2008). \emph{A multiple testing correction method for genetic association studies using correlated single nucleotide polymorphisms.} Genetic epidemiology, 32(4), 361-369.

   Galwey, N. W. (2009). \emph{A new measure of the effective number of tests, a practical tool for comparing families of non-independent significance tests.} Genetic epidemiology, 33(7), 559-568.

   Li, J., & Ji, L. (2005). \emph{Adjusting multiple testing in multilocus analyses using the eigenvalues of a correlation matrix.} Heredity, 95(3), 221-227.

   Nyholt, D. R. (2004). \emph{A simple correction for multiple testing for single-nucleotide polymorphisms in linkage disequilibrium with each other.} The American Journal of Human Genetics, 74(4), 765-769.
}
\examples{
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
meff(R, method = "nyholt")
}
\keyword{htest}
