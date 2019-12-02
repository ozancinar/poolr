\name{fisher}
\alias{fisher}
\title{Fisher's Method}
\description{Function to carry out Fisher's method.}
\usage{
fisher(p, adjust = "none", m, R, size = 10000, threshold, side = 2, batchsize, \dots)
}
\arguments{
   \item{p}{vector of length \eqn{k} with the (one- or two-sided) p-values to be combined.}
   \item{adjust}{character string to specify an adjustment method to account for dependence. Default is \code{"none"}, in which case no adjustment method is applied. Methods \code{"nyholt"}, \code{"liji"}, \code{"gao"}, or \code{"galwey"} are adjustments based on the effective number of tests (see \code{\link{meff}}). Adjustment method \code{"empirical"} uses an empirically-derived null distribution using pseudo replicates. Finally, method \code{"generalized"} uses a generalization of Fisher's method based on multivariate theory. See \sQuote{Details}.}
   \item{R}{a symmetric matrix that reflects the dependence structure among the tests. Must be specified if \code{adjust} is set to something other than \code{"none"} (or the effective number of tests is directly provided via the \code{m} argument).}
   \item{m}{optional scalar (between 1 and \eqn{k}) to manually specify the effective number of tests (instead of estimating it via one of the methods described above).}
   \item{size}{size of the empirically-derived null distribution. Can also be a numeric vector of sizes, in which case a stepwise algorithm is used. This (and the following arguments) are only relevant when \code{adjust = "empirical"}.}
   \item{threshold}{numeric vector to specify the significance thresholds for the stepwise algorithm (only relevant when \code{size} is a vector).}
   \item{side}{scalar to specify the sidedness of the \eqn{p}-values that are used to generate the empirically-derived null distribution (2, by default, for two-sided tests, 1 for one-sided tests).}
   \item{batchsize}{logical that defines how the empirically-derived null distribution is constructed. When \code{FALSE} (the default), this is done in a typically faster but more memory intensive vectorized manner. When \code{TRUE}, a typically slower but less memory intensive loop is used.}
   \item{\dots}{other arguments.}
}
\details{

   \bold{Fisher's Method}

   By default (i.e., when \code{adjust = "none"}), the function applies Fisher's (1932) method to the \eqn{p}-values passed to the function. The test statistic is then computed with \deqn{X^2 = -2 \sum_{i = 1}^k \ln(p_i)} where \eqn{p_1, p_2, \ldots, p_k} are the individual \eqn{p}-values of the \eqn{k} hypothesis tests. Under the joint null hypothesis, the test statistic follows a chi-squared distribution with \eqn{2k} degrees of freedom which is used to compute the combined \eqn{p}-value.

   Fisher's method assumes that the \eqn{p}-values to be combined are independent. If this is not the case, the method can either be conservative (not reject often enough) or liberal (reject too often), depending on the correlation structure among the tests. In this case, one can adjust Fisher's method to account for such dependence.

   \bold{Adjustment Based on the Effective Number of Tests}

   The method can be adjusted based on an estimate of the effective number of tests. For this, set \code{adjust} to \code{"nyholt"}, \code{"liji"}, \code{"gao"} or \code{"galwey"} (see \code{\link{meff}} for details on these methods for estimating the effective number of tests) and argument \code{R} to a correlation matrix that reflects the dependence structure among the tests.

   There is no general solution for constructing such a correlation matrix, as this depends on the type of test that generated the \eqn{p}-values and the sidedness of these tests. In case that the \eqn{p}-values are obtained from tests whose test statistics can be assumed to follow (at least asymptotically) a multivariate normal distribution and a matrix is available that reflects the correlations among the test statistics, then the \code{\link{mvnconv}} function can be used to convert this correlation matrix into the correlations among the (one- or two-sided) \eqn{p}-values, which can then be passed to the \code{R} argument. See \sQuote{Examples}.

   Once the effective number of tests, \eqn{m}, is estimated based on \code{R} using one of the four methods described above, the test statistic of Fisher's method can be modified with \deqn{\tilde{X}^2 = \frac{m}{k} \times X^2} which is then assumed to follow a chi-squared distribution with \eqn{2m} degrees of freedom.

   Alternatively, one can also directly specify the effective number of tests via the \code{m} argument (e.g., if some other method not implemented in the \pkg{poolr} package is used to estimate the effective number of tests). In that case, argument \code{R} is irrelevant and doesn't need to be specified.

   \bold{Adjustment Based on an Empirically-Derived Null Distribution}

   When \code{adjust="empirical"}, the combined \eqn{p}-value is computed based on an empirically-derived null distribution of Fisher's method using pseudo replicates. This is appropriate if the test statistics that generated the \eqn{p}-values to be combined can be assumed to follow (at least asymptotically) a multivariate normal distribution and a matrix is available that reflects the correlations among the test statistics (which is specified via the \code{R} argument). In this case, test statistics are repeatedly simulated from a multivariate normal distribution under the joint null hypothesis, converted into one- or two-sided \eqn{p}-values (depending on the \code{side} argument), and Fisher's method is applied. Repeating this process \code{size} times yields a null distribution based on which the combined \eqn{p}-value can be computed, or more precisely, estimated, since repeated applications of this method will yield (slightly) different results.

   To obtain a stable estimate of the combined \eqn{p}-value, especially in case the combined \eqn{p}-value ends up being quite low (i.e., close to 0), \code{size} should be set to a large value (the default is \code{10000}, but this might need to be increased for a more precise estimate).

   Generating the null distribution can be done in two ways. When \code{emploop=FALSE} (the default), this is done in a typically faster but more memory intensive vectorized manner, which requires generating a matrix with dimensions \code{[size,k]}. For large \code{size} and/or \code{k}, allocating the memory for this matrix might not be possible. In that case, one can set \code{emploop=TRUE}, which uses a typically slower but less memory intensive loop.

   One can also specify a vector for the \code{size} argument, in which case one must also specify a corresponding vector for the \code{threshold} argument. In that case, a stepwise algorithm is used that proceeds as follows (see also Liu et al., 2010). For \code{j = 1, ..., length(size)},
   \enumerate{
      \item estimate the combined \eqn{p}-value based on \code{size[j]}
      \item if the combined \eqn{p}-value is \eqn{\ge} than \code{threshold[j]}, stop (and report the combined \eqn{p}-value), otherwise go back to 1.
   }
   By setting \code{size} to increasing values (e.g., \code{size = c(1000, 100000, 10000000)} and \code{threshold} to decreasing values (e.g., \code{threshold = c(.10, .01, 0)}), one can quickly obtain a fairly accurate estimate of the combined \eqn{p}-value if it is far from significant (e.g., \eqn{\ge .10}), but hone in on a more accurate estimate for a combined \eqn{p}-value that is closer to 0. Note that the last value of \code{threshold} should, and is forced to be, 0 so that the algorithm is guaranteed to terminate (hence, one can also leave out the last value of \code{threshold}, so \code{threshold = c(.10, .01)} would also work in the example above).

   \bold{Adjustment Based on Multivariate Theory}

   Finally, the test-statistics can be calculated under dependence which provides the combined \eqn{p}-value that accounts for the correlation structure (Brown, 1975). This requires the calculation of the correlations among the test-statistics which can be done by using \code{mvnconv()}. For details, see \code{?mvnconv)} and examples below.
}
\value{
   The function returns the (pooled) p-value based on Fisher's method and the selected adjustment method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Brown, MB (1975): \emph{400: A method for combining non-independent, one-sided tests of significance}. Biometrics, 987-992.

   Fisher, R. A. (1932). \emph{Statistical Methods for Research Workers} (4th ed.). Edinburgh: Oliver and Boyd.

   Liu, J. Z., McRae, A. F., Nyholt, D. R., Medland, S. E., Wray, N. R., Brown, K. M., AMFS Investigators, Hayward, N. K., Montgomery, G. W., Visscher, P. M., Martin, N. G., & Macgregor, S. (2010). A versatile gene-based test for genome-wide association studies. American Journal of Human Genetics, 87(1), 139-145.
}
\examples{
# create a pseudo-correlation matrix for 10 test-statistics
r <- matrix(0.5, 10, 10)
diag(r) <- 1

# generate multivariate normal test statistics with given correlation structure
set.seed(123)
ts <- mvrnorm(n = 1, mu = rep(0, 10), Sigma = r)

# compute one-sided p-values from the test statistics generated above
p1 <- pnorm(ts)

# compute two-sided p-values from the test-statistics generated above
p2 <- 2 * pnorm(abs(ts), lower.tail = FALSE)

# apply standard Fisher's method
fisher(p1)
fisher(p2)

# adjustment using a PCA-based method (Li & Ji)
fisher(p1, adjust = "liji", R = mvnconv(r, side = 1, target = "p", cov2cor = TRUE))
fisher(p2, adjust = "liji", R = mvnconv(r, side = 2, target = "p", cov2cor = TRUE))

# use empirically-derived null distributions
fisher(p1, adjust = "empirical", R = r, side = 1)
fisher(p2, adjust = "empirical", R = r) # side = 2 is default

# using empirically-derived null distributions with a stepwise algorithm.
fisher(p2 / 5, adjust = "empirical", R = r,
       size = c(1000, 10000, 100000), threshold = c(0.5, 0.5), verbose = TRUE)

# compute combined p-value under dependence
fisher(p1, adjust = "generalized", R = mvnconv(r, side = 1, target = "m2lp"))
fisher(p2, adjust = "generalized", R = mvnconv(r, side = 2, target = "m2lp"))
}
\keyword{htest}
