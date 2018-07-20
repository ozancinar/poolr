\name{bonferroni}
\alias{bonferroni}
\title{Bonferroni Correction}
\description{Function to combine the p-values by using the Bonferroni correction.}
\usage{
bonferroni(p, adjust = "none", m, R,
           size = 10000, seed, type = 2, emp.loop = FALSE, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. As default, the Bonferroni correction is applied to the minimum p-value. Optionally, number of significant tests can be used instead of the Bonferroni correction by setting this argument to one of the four PCA-based methods, "nyholt", "liji", "gao" or "galwey". Finally, this argument can be set to "empirical" to use the empircally-derived distributions under the null hypothesis.}
   \item{m}{an optional argument that can be used to set the effective number of tests to a user-defined value. This argument can only be used when the "adjust" argument is set to one of the four PCA-based methods.}
   \item{R}{a symmetric matrix that includes pairwise correlations between the p-values.}
   \item{size}{size of the emprical distribution if the empirically-derived null hypothesis is going to be used.}
   \item{seed}{seed for the random number generation if the empirically-derived null hypothesis is going to be used.}
   \item{type}{type of the p-values used in the tests. As default, it is set to 2 for two-sided tests.}
   \item{emp.loop}{a logical argument to decide to generate the empirical distribution in a loop. It is set to FALSE as default. However, if the size of the empirical distribution to be generated is too large, we recommend you to switch emp.loop = TRUE.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to combine a set of p-values by using the Bonferroni correction and the strongest evidence, i.e. the smallest p-value.
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
}
\examples{
p <- c(.013, .022, .124, .278)
R <- matrix(c(1, 0.8, 0.5, -0.3,
              0.8, 1, 0.2, 0.4,
              0.5, 0.2, 1, -0.7,
              -0.3, 0.4, -0.7, 1),
              nrow = 4, byrow = TRUE)
bonferroni(p)
bonferroni(p, adjust = "nyholt", R = R)
bonferroni(p, adjust = "empirical", R = R)
}
\keyword{htest}
