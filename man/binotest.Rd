\name{binotest}
\alias{binotest}
\title{Binomial Test}
\description{Function to carry out the binomial test to combine p-values.}
\usage{
binotest(p, adjust = "none", m, R, alpha = 0.05,
         size = 10000, seed, type = 2, emp.loop = FALSE, emp.step, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. Can be "none", as default, to use the default binomial test. Or the test can be adjusted with PCA-based methods by calculating the number of effective tests. One of the four PCA-based methods can be used by setting this argument to one of "nyholt", "liji", "gao" or "galwey". Finally, this argument can be set to "empirical" to use empirically-derived distributions under the null hypothesis.}
   \item{m}{an optional argument that can be used to set the effective number of tests to a user-defined value. This argument can only be used when the "adjust" argument is set to one of the four PCA-based methods.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{alpha}{significance level of the test. It is set to 0.05 as default.}
   \item{size}{size of the emprical distribution if the empirically-derived null hypothesis is going to be used.}
   \item{seed}{seed for the random number generation if the empirically-derived null hypothesis is going to be used.}
   \item{type}{type of the p-values used in the tests. As default, it is set to 2 for two-sided tests.}
   \item{emp.loop}{a logical argument to decide to generate the empirical distribution in a loop. It is set to FALSE as default. However, if the size of the empirical distribution to be generated is too large, we recommend you to switch emp.loop = TRUE.}
   \item{emp.step}{an argument to apply a stepwise algortihm for empirically-derived null distributions. See details for further information.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to use binomial test to combine p-values. The user can choose to use default parameters for the binomial test depending on the number of p-values to be combined by leaving "adjust" argument as default. Or, the user can choose to adjust the 'n' parameter of the binomial distribution with PCA-based methods to consider the dependencies between the tests.
   
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
p <- c(.013, .022, .124, .278)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
binotest(p)
binotest(p, adjust = "liji", R = R)
binotest(p, adjust = "empirical", R = R)
}
\keyword{htest}
