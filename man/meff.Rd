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
   The function can be used to calculate the effective number of tests based on one of the four methods. \code{poolR} provides four PCA-based methods to calculate the effective number of tests. In order to find the effective number of tests with any of these methods, first the eigen values based on the correlation matrix are computed where \eqn{\lambda_i} denotes the \eqn{i}th eigen-value. Then, the computations of effective number of tests based on the eigen-values are shown below. 

   - Nyholt (2004): \eqn{m_\textrm{eff} = 1 + (k - 1) \left(1 - \frac{\mbox{Var}(\lambda)}{k}\right)}

   where \eqn{\mbox{Var}(\lambda)} is the variance of the \eqn{k} eigenvalues.

   - Li & Ji (2005): \eqn{m_\textrm{eff} = \sum_{i = 1}^k f(|\lambda_i|)}

   where \eqn{f(x) = I(x \geq 1) + (x - \lfloor x \rfloor)} and \eqn{\lfloor.\rfloor} is the floor function.

   - Gao et al. (2008): \eqn{m_\textrm{eff} = \min(x) \; \mbox{such that} \; \frac{\sum_{i = 1}^x \lambda_{(i)}}{\sum_{i = 1}^k \lambda_{(i)}} > C}

   where \eqn{\lambda_{(1)}, \lambda_{(2)}, \ldots, \lambda_{(k)}} denote the eigen values in a decreasing order and \eqn{C} is a pre-defined parameter which is set to 0.995. 

   - Galwey (2009): \eqn{m_\textrm{eff} = \frac{\left(\sum_{i = 1}^k \sqrt{\lambda_i}\right)^2}{\sum_{i = 1}^k \lambda_i}}

   In real life studies, researchers probably have the correlations among the interchanging elements or the test-statistics. For techniques to estimate the correlations among the \eqn{p}-values based on the correlation matrix at hand, please see \code{?mvnconv}. 

}
\value{
   The function returns the effective number of tests.
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
