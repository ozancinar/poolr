\name{bonferroni}
\alias{bonferroni}
\title{Bonferroni Correction}
\description{Function to combine the \eqn{p}-values by using the Bonferroni correction.}
\usage{
bonferroni(p, adjust = "none", m, R,
           size = 10000, seed, side = 2, emp.loop = FALSE, emp.step, \dots)
}
\arguments{
   \item{p}{vector of the (one- or two-sided) \eqn{p}-values to be combined.}
   \item{adjust}{adjustment method while combining the \eqn{p}-values. It is set to \code{"none"} as default which assumes the p-values are independent. The process of combining \eqn{p}-values can be modified with PCA-based methods by setting this argument to one of the four such methods which are \code{"nyholt"}, \code{"liji"}, \code{"gao"}, or \code{"galwey"} (for details of PCA-based methods, see \code{?meff}). Or, this argument can be set to \code{"empirical"} to use empirically-derived distributions under the null hypothesis.}
   \item{m}{an optional argument that can be used to set the effective number of tests to a user-defined value. Once this argument is set to a value, the function will not use \code{adjust} and \code{R} arguments.}
   \item{R}{a symmetric matrix that stores pairwise correlations among the \eqn{p}-values. The user must submit a matrix when the \code{adjust} set to a value other than \code{"none"} or an effective number of tests is provided with \code{m} argument. To estimate the correlations among the \eqn{p}-values for different methods, please see \code{?mvnconv}.}
   \item{size}{size of the emprical distribution if the empirically-derived null hypothesis is going to be used.}
   \item{seed}{seed for the random number generation if the empirically-derived null hypothesis is going to be used.}
   \item{side}{type of the \eqn{p}-values that are generated in empirical distribution. As default, it is set to 2 for two-sided tests. For one-sided \eqn{p}-values it should be set to 1.}
   \item{emp.loop}{a logical argument to decide to generate the empirical distribution in a loop. It is set to \code{FALSE} as default. However, if the size of the empirical distribution to be generated is so large that it leads to memory problem, we recommend you to switch \code{emp.loop = TRUE}.}
   \item{emp.step}{an argument to apply a stepwise algortihm for empirically-derived null distributions. See details for further information.}
   \item{\dots}{other arguments.}
}
\details{
  The function can be used to combine a set of \eqn{p}-values by utilizing the Bonferroni correction which uses the strongest evidence, i.e., the smallest \eqn{p}-value among individual \eqn{p}-values denoted by \eqn{p_1, p_2, \ldots, p_k} where \eqn{k} is the number of individual \eqn{p}-values. This is done as follows under complete independence:

  \eqn{p = \min(1, \min(p_1, \ldots, p_k) \times k)}

  where \eqn{p} is the combined \eqn{p}-value.

  When the hypothesis tests from which the individual \eqn{p}-values have been obtained are not completely independent, this process can be modified with several techniques. First, users can utilize the effective number of tests by setting \code{adjust} argument to one of the four PCA-based methods \code{"nyholt"}, \code{"liji"}, \code{"gao"} or \code{"galwey"} (see \code{?meff} for details) and providing the correlation among the \eqn{p}-values (see \code{?mvnconv} for more details). Alternatively, users can provide an effective number of tests which is computed by another method by submitting this value to \code{m} argument. Note that once a user-defined effective number of tests is provided with \code{m} argument, he/she does not have to provide a correlation matrix to \code{R} argument anymore. Once the effective number of tests, \eqn{m_\textrm{eff}}, is computed, the combined \eqn{p}-value can be modified with it as:

  \eqn{p = \min(1, \min(p_1, \ldots, p_k) \times m_\textrm{eff}}.

  Alternatively, researchers can utilize the empirically-derived null distributions in order to calculate the combined \eqn{p}-value under a given correlation structure. This can be done by setting \code{adjust} to \code{"empirical"} and providing the correlation matrix to \code{R} argument. Then, the function will generate an empirical distribution of the combined \eqn{p}-values under the null hypothesis and compute a combined \eqn{p}-value using the percentiles of this distribution. In order to quicken the computations, researchers can use a stepwise algorithm. For details, please see below.

  \bold{Stepwise Algorithm for Empirically-Derived Null Distributions}

   As an option the \code{empirical} adjustment can be applied in a stepwise manner. With this approach, the function expects the user to set different number of sample size values for the empirical data set, and respective significance thresholds for these sample size values except for the largest sample size. This information should be given to the function as a list, containing two vectors for sample size and significance threshold values. Optimally, these vectors should be named \code{size} for sample size values and \code{thres} for significance threshold values, where the length of \code{thres} should be one less number than the length of \code{size} (although these names are going to be set by the function if they are not given by the user, as long as the list object has appropriate vectors). Then, the function will follow a loop over the sample size values in the order that is given in the \code{size} vector except for the last sample size value. Note that the algorithm will follow the order that is given by the user. Therefore, we advise you to give the sample size values in the \code{size} vector in an ascending order and their respective significance threshold values to benefit this algorithm optimally. Until the last sample size value, the function will generate an empirical distribution with the given sample size value in the order. This empirical distribution is going to be used to calculate the combined p-value. If the combined p-value is larger than the respective significance threshold, the algorithm will stop and give this p-value and the respective test statistic as the result. If the combined p-value is smaller than the respective significance threshold, then the algorithm will continue with the next sample size value until the second last sample size value in the set. If none of the combined p-values from these sample size values are below the respective significance threshold values, then the algorithm will generate one last empirical distribution with the last sample size value in the list and calculate the combined p-value and the test statistic based on this empirical distribution. For another explanation please see Liu et al. (2010).
}
\value{
   The function returns the (pooled) p-value based on the Bonferroni correction.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Bland, J. M., & Altman, D. G. (1995). \emph{Multiple significance tests: The Bonferroni method}. British Medical Journal, 310(6973), 170.

   Jimmy Z. Liu, Allan F. Mcrae, Dale R. Nyholt, Sarah E. Medland, Naomi R. Wray, Kevin M. Brown, AMFS Investigators, Nicholas K. Hayward, Grant W. Montgomery, Peter M. Visscher, Nicholas G. Martin, and Stuart Macgregor (2010). \emph{A versatile gene-based test for genome-wide association studies}. The American Journal of Human Genetics, 87(1), 139-145.
}
\examples{
# create a pseudo-correlation matrix for 10 test-statistics
r <- matrix(0.3, 10, 10)
diag(r) <- 1

# generating multivariate normal test-statistics with the given correlation structure
set.seed(123)
ts <- mvrnorm(n = 1, mu = rep(0, 10), Sigma = r)

# computing individual one-sided p-values from the test-statistics generated above
p1 <- pnorm(ts)

# computing individual two-sided p-values from the test-statistics generated above
p2 <- 2 * pnorm(abs(ts), lower.tail = FALSE)

# default Bonferroni method
bonferroni(p1)
bonferroni(p2)

# using a PCA-based method (Li & Ji)
bonferroni(p1, adjust = "liji", R = mvnconv(R = r, side = 1, target = "p"))
bonferroni(p2, adjust = "liji", R = mvnconv(R = r, side = 2, target = "p"))

# using empirically-derived null distributions
bonferroni(p1, adjust = "empirical", R = r, side = 1)
bonferroni(p2, adjust = "empirical", R = r)  # side = 2 as default.
bonferroni(p2, adjust = "empirical", R = r, seed = 321)  # setting the seed for
                                                             # reproducibility.
bonferroni(p2, adjust = "empirical", R = r,
           emp.step = list(size = c(1000, 10000, 100000), thres = c(0.5, 0.05)))
}
\keyword{htest}
