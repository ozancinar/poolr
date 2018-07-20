\name{binotest}
\alias{binotest}
\title{Binomial Test}
\description{Function to carry out the binomial test to combine p-values.}
\usage{
binotest(p, adjust = "none", m, R, alpha = 0.05,
         size = 10000, seed, type = 2, emp.loop = FALSE, \dots)
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
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to use binomial test to combine p-values. The user can choose to use default parameters for the binomial test depending on the number of p-values to be combined by leaving "adjust" argument as default. Or, the user can choose to adjust the 'n' parameter of the binomial distribution with PCA-based methods to consider the dependencies between the tests.
}
\value{
   The function returns the (pooled) p-value calculated with binomial distribution and the selected adjustment method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
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
