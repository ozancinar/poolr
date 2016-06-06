\name{meff}
\alias{meff}
\title{Number of Effective Tests}
\description{Function to calculate the number of effective tests based on one of the four methods.}
\usage{
meff(R, method, \dots)
}
\arguments{
   \item{R}{matrix that include pairwise correlations between the p-values.}
   \item{method}{method to be used to calculate the number of effective tests. Can be one of "nyholt", "li.ji", "gao", or "galwey".}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to calculate the number of effective tests based on one of the four methods.
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
meff(R = R, method = "nyholt")
}
\keyword{htest}
