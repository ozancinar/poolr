\name{binotest}
\alias{binotest}
\title{Binomial Test}
\description{Function to combine \eqn{p}-values using the binomial distribution.}
\usage{
binotest(p, adjust = "none", m, R, size = 10000, threshold, side = 2, batchsize, \dots)
}
\arguments{
   \item{p}{vector of the (one- or two-sided) \eqn{p}-values to be combined.}
   \item{adjust}{adjustment method while combining the \eqn{p}-values. It is set to \code{"none"} as default which assumes the p-values are independent. The process of combining \eqn{p}-values can be modified with PCA-based methods by setting this argument to one of the four such methods which are \code{"nyholt"}, \code{"liji"}, \code{"gao"}, or \code{"galwey"} (for details of PCA-based methods, see \code{?meff}). Or, this argument can be set to \code{"empirical"} to use empirically-derived distributions under the null hypothesis.}
   \item{m}{an optional argument that can be used to set the effective number of tests to a user-defined value. Once this argument is set to a value, the function will not use \code{adjust} and \code{R} arguments.}
   \item{R}{a symmetric matrix that stores pairwise correlations among the \eqn{p}-values. The user must submit a matrix when the \code{adjust} set to a value other than \code{"none"} or an effective number of tests is provided with \code{m} argument. To estimate the correlations among the \eqn{p}-values for different methods, please see \code{?mvnconv}.}
   \item{size}{size of the emprical distribution if the empirically-derived null hypothesis is going to be used. For a stepwise algorithm this qrgument can be given as a numeric vector to specify the sample sizes at different steps.}
   \item{threshold}{a numeric vector to specify significance thresholds for the stepwise algorithm. This option is used when \code{size} is set to a vector. Then, \code{threshold} should be set to a vector with a length of one less than the length of \code{size}. If these two vectors have the same lengths, the last threshold in \code{threshold} is going to be ignored. Alternatively, \code{threshold} can be a single value which will be used for the signficance threshold for all the steps except for the last step.}
   \item{side}{type of the \eqn{p}-values that are generated in empirical distribution. As default, it is set to 2 for two-sided tests. For one-sided \eqn{p}-values it should be set to 1.}
   \item{batchsize}{a logical argument to decide to generate the empirical distribution in a loop. It is set to \code{FALSE} as default. However, if the size of the empirical distribution to be generated is so large that it leads to memory problem, we recommend you to switch \code{emploop = TRUE}.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to combine p-values based on a binomial distribution. It is known that \eqn{p}-values follow Uniform(0, 1) under the null hypothesis. Therefore, rejecting the null hypothesis (i.e., a Type I Error) is equal to \eqn{\alpha}. A combined \eqn{p}-value, \eqn{p}, that indicates an ``excess significance'' can then be calculated under independence as:

   \eqn{p = \sum_{x = r}^m {m \choose x} \alpha^x (1 - \alpha)^{m - x}}

   where \eqn{r} is the number of individual \eqn{p}-values that are smaller than \eqn{\alpha}.

   When the hypothesis tests from which the individual \eqn{p}-values have been obtained are not completely independent, this process can be modified with several techniques. First, users can utilize the effective number of tests by setting \code{adjust} argument to one of the four PCA-based methods \code{"nyholt"}, \code{"liji"}, \code{"gao"} or \code{"galwey"} (see \code{?meff} for details) and providing the correlation among the \eqn{p}-values (see \code{?mvnconv} for more details). Alternatively, users can provide an effective number of tests which is computed by another method by submitting this value to \code{m} argument. Note that once a user-defined effective number of tests is provided with \code{m} argument, he/she does not have to provide a correlation matrix to \code{R} argument anymore. Once the effective number of tests, \eqn{m_\textrm{eff}}, is computed, the combined \eqn{p}-value can be modified with it as:

   \eqn{p = \sum_{x=\tilde{r}}^{m_{\textrm{eff}}} {m_{\textrm{eff}} \choose x} \alpha_p^x (1-\alpha_p)^{m_{\textrm{eff}}-x}}

   where \eqn{\tilde{r} = \lfloor r \frac{k_{\textrm{eff}}}{k} \rfloor}.

   Alternatively, researchers can utilize the empirically-derived null distributions in order to calculate the combined \eqn{p}-value under a given correlation structure. This can be done by setting \code{adjust} to \code{"empirical"} and providing the correlation matrix to \code{R} argument. Then, the function will generate an empirical distribution of the combined \eqn{p}-values under the null hypothesis and compute a combined \eqn{p}-value using the percentiles of this distribution. In order to quicken the computations, researchers can use a stepwise algorithm. For details, please see below.

   \bold{Stepwise Algorithm for Empirically-Derived Null Distributions}

   As an option the \code{empirical} adjustment can be applied in a stepwise manner. With this approach, the function expects the user to set different number of sample size values for the empirical data set, and respective significance thresholds for these sample size values except for the largest sample size. This information should be given to the function as a list, containing two vectors for sample size and significance threshold values. Optimally, these vectors should be named \code{size} for sample size values and \code{thres} for significance threshold values, where the length of \code{thres} should be one less number than the length of \code{size} (although these names are going to be set by the function if they are not given by the user, as long as the list object has appropriate vectors). Then, the function will follow a loop over the sample size values in the order that is given in the \code{size} vector except for the last sample size value. Note that the algorithm will follow the order that is given by the user. Therefore, we advise you to give the sample size values in the \code{size} vector in an ascending order and their respective significance threshold values to benefit this algorithm optimally. Until the last sample size value, the function will generate an empirical distribution with the given sample size value in the order. This empirical distribution is going to be used to calculate the combined p-value. If the combined p-value is larger than the respective significance threshold, the algorithm will stop and give this p-value and the respective test statistic as the result. If the combined p-value is smaller than the respective significance threshold, then the algorithm will continue with the next sample size value until the second last sample size value in the set. If none of the combined p-values from these sample size values are below the respective significance threshold values, then the algorithm will generate one last empirical distribution with the last sample size value in the list and calculate the combined p-value and the test statistic based on this empirical distribution. For another explanation please see Liu et al. (2010).
}
\value{
   The function returns the (pooled) p-value calculated with binomial distribution and the selected adjustment method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Jimmy Z. Liu, Allan F. Mcrae, Dale R. Nyholt, Sarah E. Medland, Naomi R. Wray, Kevin M. Brown, AMFS Investigators, Nicholas K. Hayward, Grant W. Montgomery, Peter M. Visscher, Nicholas G. Martin, and Stuart Macgregor (2010). \emph{A versatile gene-based test for genome-wide association studies}. The American Journal of Human Genetics, 87(1), 139-145.

   Wilkinson, B. (1951). \emph{A statistical consideration in psychological research}. Psychological Bulletin, 48(2), 156-158.
}

\examples{
# create a pseudo-correlation matrix for 10 test-statistics.
r <- matrix(0.3, 10, 10)
diag(r) <- 1

# generating multivariate normal test-statistics with the given
# correlation structure.
set.seed(123)
ts <- mvrnorm(n = 1, mu = rep(0, 10), Sigma = r)

# computing individual one-sided p-values from the test-statistics
# generated above.
p1 <- pnorm(ts)

# computing individual two-sided p-values from the test-statistics
# generated above.
p2 <- 2 * pnorm(abs(ts), lower.tail = FALSE)

# default binomial test method.
binotest(p = p1)
binotest(p = p2)

# using a PCA-based method (Li & Ji).
binotest(p = p1, adjust = "liji", R = mvnconv(R = r, side = 1, target = "p"))
binotest(p = p2, adjust = "liji", R = mvnconv(R = r, side = 2, target = "p"))

# using empirically-derived null distributions.
binotest(p = p1, adjust = "empirical", R = r, side = 1)
binotest(p = p2, adjust = "empirical", R = r)  # side = 2 as default.

# using empirically-derived null distributions with a stepwise algorithm.
binotest(p = p2 / 5, adjust = "empirical", R = r,
         size = c(1000, 10000, 100000), thres = c(0.5, 0.05), verbose = TRUE)
}
\keyword{htest}
