\name{empirical}
\alias{empirical}
\title{Simulate Empirically-Derived Null Distributions}
\description{Function to simulate empirically-derived null distributions of various methods for combining \mjseqn{p}-values using pseudo replicates.\loadmathjax}
\usage{
empirical(R, method, side = 2, size = 10000, batchsize, \dots)
}
\arguments{
   \item{R}{a \mjeqn{k \times k}{k * k} symmetric matrix that contains the correlations among the test statistics.}
   \item{method}{character string to specify for which method to simulate the null distribution (either \code{"fisher"}, \code{"stouffer"}, \code{"invchisq"}, \code{"binotest"}, \code{"bonferroni"}, or \code{"tippett"}).}
   \item{side}{scalar to specify the sidedness of the \mjseqn{p}-values that are used to simulate the null distribution (2, by default, for two-sided tests; 1 for one-sided tests).}
   \item{size}{size of the empirically-derived null distribution that should be generated.}
   \item{batchsize}{optional scalar to specify the batch size for generating the null distribution. When unspecified (the default), this is done in a single batch.}
   \item{\dots}{other arguments.}
}
\details{

   This function simulates the null distribution of a particular method for combining \mjseqn{p}-values when the test statistics that generate the \mjseqn{p}-values to be combined can be assumed to follow a multivariate normal distribution and a matrix is available that reflects the correlations among the test statistics (which is specified via the \code{R} argument). In this case, test statistics are repeatedly simulated from a multivariate normal distribution under the joint null hypothesis, converted into one- or two-sided \mjseqn{p}-values (depending on the \code{side} argument), and the chosen method is applied. Repeating this process \code{size} times yields the null distribution.

   If \code{batchsize} is unspecified, the null distribution is simulated in a single batch, which requires temporarily storing a matrix with dimensions \code{[size,k]}. When \code{size} and/or \code{k} is large, allocating the memory for this matrix might not be possible. Instead, one can specify a \code{batchsize} value, in which case a matrix with dimensions \code{[batchsize,k]} is repeatedly simulated until the desired size of the null distribution has been obtained.

}
\value{
   A vector of combined \mjseqn{p}-values as simulated under the joint null hypothesis for a given method.
}
\note{
   The \code{R} matrix must be positive definite. If it is not, the function uses \code{\link[Matrix]{nearPD}} to find the nearest positive definite matrix (Higham, 2002) before simulating the null distribution.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Higham, N. J. (2002). Computing the nearest correlation matrix: A problem from finance. \emph{IMA Journal of Numerical Analysis, 22}(3), 329-343.
}
\examples{
# create an example correlation matrix with constant positive correlations
R <- matrix(0.6, nrow = 10, ncol = 10)
diag(R) <- 1

# generate null distribution for Fisher's method (setting the seed for reproducibility)
set.seed(1234)
psim <- empirical(R, method = "fisher")

# Fisher's method is liberal in this scenario (i.e., its actual Type I error
# rate is around .14 instead of the nominal significance level of .05)
mean(psim <= .05)

# estimate the actual Type I error rate of the other methods in this scenario
psim <- empirical(R, method = "stouffer")
mean(psim <= .05)
psim <- empirical(R, method = "invchisq")
mean(psim <= .05)
psim <- empirical(R, method = "binotest")
mean(psim <= .05)
psim <- empirical(R, method = "bonferroni")
mean(psim <= .05)
psim <- empirical(R, method = "tippett")
mean(psim <= .05)

# Stouffer's and the inverse chi-square method also have clearly inflated
# Type I error rates and the binomial test just barely. As expected, the
# Bonferroni method is overly conservative and so is Tippett's method.
}
\keyword{htest}
