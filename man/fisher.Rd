\name{fisher}
\alias{fisher}
\title{Fisher's Method}
\description{Function to carry out Fisher's method.}
\usage{
fisher(p, adjust = "none", m, R,
       size = 10000, seed, type = 2, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. Can be "none", as default, to use the default Fisher's method. Furthermore, Fisher's method can be adjusted with a various ways to consider the dependencies between the individual p-values. First, this argument can be set to one of "nyholt", "liji", "gao", or "galwey" to use the PCA-based methods by calculating the number of effective tests. Second, it can be set to "empirical" to use empirically-derived distributions under the null. Finally, this argument can be set to "brown" to use the Brown's method.}
   \item{m}{an optional argument that can be used to set the effective number of tests to a user-defined value. This argument can only be used when the "adjust" argument is set to one of the four PCA-based methods.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{size}{size of the emprical distribution if the empirically-derived null hypothesis is going to be used.}
   \item{seed}{seed for the random number generation if the empirically-derived null hypothesis is going to be used.}
   \item{type}{type of the p-values used in the tests. As default, it is set to 2 for two-sided tests.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to apply Fisher's method. The user can choose to use default Fisher's method by leaving "adjust" argument as default. Or, the user can choose to adjust the Fisher's method with PCA-based methods to consider the dependencies between the tests.
}
\value{
   The function returns the (pooled) p-value based on Fisher's method and the selected adjustment method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Fisher, R. A. (1932). \emph{Statistical methods for research workers} (4th ed.). Oliver and Boyd.

   Brown, MB (1975): \emph{400: A method for combining non-independent, one-sided tests of significance}. Biometrics, 987-992.
}
\examples{
p <- c(.013, .022, .124, .278)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
fisher(p)
fisher(p, adjust = "gao", R = R)
fisher(p, adjust = "empirical", R = R)
fisher(p, adjust = "brown", R = R)
}
\keyword{htest}
