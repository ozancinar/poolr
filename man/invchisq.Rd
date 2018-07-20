\name{invchisq}
\alias{invchisq}
\title{Inverse Chi-Squared Method}
\description{Function to carry out inverse Chi-Squared method.}
\usage{
invchisq(p, adjust = "none", m, R,
         size = 10000, seed, type = 2, emp.loop = FALSE, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. Can be "none", as default, to use the default inverse Chi-Squared method. This method can be adjusted with a various ways to consider the dependencies between the individual p-values. First, this argument can be set to one of "nyholt", "liji", "gao", or "galwey" to use the PCA-based methods by calculating the number of effective tests. Or, it can be set to "empirical" to use an empirically-derived distribution of the test statistic under the null hypothesis.}
   \item{m}{an optional argument that can be used to set the effective number of tests to a user-defined value. This argument can only be used when the "adjust" argument is set to one of the four PCA-based methods.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{size}{size of the emprical distribution if the empirically-derived null hypothesis is going to be used.}
   \item{seed}{seed for the random number generation if the empirically-derived null hypothesis is going to be used.}
   \item{type}{type of the p-values used in the tests. As default, it is set to 2 for two-sided tests.}
   \item{emp.loop}{a logical argument to decide to generate the empirical distribution in a loop. It is set to FALSE as default. However, if the size of the empirical distribution to be generated is too large, we recommend you to switch emp.loop = TRUE.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to apply inverse Chi-Squared method. The user can choose to use default inverse Chi-Squared method by leaving "adjust" argument as default. The user can choose to adjust the method with PCA-based methods to consider the dependencies between the tests. Or, the user can apply empirically-derived distribution of the test statistic under the null hypothesis to consider the correlations among the tests.
}
\value{
   The function returns the (pooled) p-value based on inverse Chi-Squared method and the selected adjustment method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Liu, J. Z., Mcrae, A. F., Nyholt, D. R., Medland, S. E., Wray, N. R., Brown, K. M., AMFS Investigators, Hayward, N. K., Montgomery, G. W., Visscher, P. M., Martin, N. G. & Macgregor, S. (2010). \emph{A versatile gene-based test for genome-wide association studies.} The American Journal of Human Genetics, 87(1), 139-145.
}
\examples{
p <- c(.013, .022, .124, .278)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
invchisq(p)
invchisq(p, adjust = "galwey", R = R)
invchisq(p, adjust = "empirical", R = R)
}
\keyword{htest}

