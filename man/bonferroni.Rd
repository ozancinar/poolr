\name{bonferroni}
\alias{bonferroni}
\title{Bonferroni Method}
\description{Function to carry out the Bonferroni method.}
\usage{
bonferroni(p, adjust = "none", R, m, size = 10000, threshold, side = 2, batchsize, \dots)
}
\arguments{
   \item{p}{vector of length \eqn{k} with the (one- or two-sided) p-values to be combined.}
   \item{adjust}{character string to specify an adjustment method to account for dependence. The default is \code{"none"}, in which case no adjustment is applied. Methods \code{"nyholt"}, \code{"liji"}, \code{"gao"}, or \code{"galwey"} are adjustments based on an estimate of the effective number of tests (see \code{\link{meff}}). Adjustment method \code{"empirical"} uses an empirically-derived null distribution using pseudo replicates. See \sQuote{Details}.}
   \item{R}{a \eqn{k \times k} symmetric matrix that reflects the dependence structure among the tests. Must be specified if \code{adjust} is set to something other than \code{"none"}. See \sQuote{Details}.}
   \item{m}{optional scalar (between 1 and \eqn{k}) to manually specify the effective number of tests (instead of estimating it via one of the methods described above).}
   \item{size}{size of the empirically-derived null distribution. Can also be a numeric vector of sizes, in which case a stepwise algorithm is used. This (and the following arguments) are only relevant when \code{adjust = "empirical"}.}
   \item{threshold}{numeric vector to specify the significance thresholds for the stepwise algorithm (only relevant when \code{size} is a vector).}
   \item{side}{scalar to specify the sidedness of the \eqn{p}-values that are used to simulate the null distribution (2, by default, for two-sided tests; 1 for one-sided tests).}
   \item{batchsize}{optional scalar to specify the batch size for generating the null distribution. When unspecified (the default), this is done in a single batch.}
   \item{\dots}{other arguments.}
}
\details{

   \bold{Bonferroni Method}

   By default (i.e., when \code{adjust = "none"}), the function applies the Bonferroni method to the \eqn{p}-values. Letting \eqn{p_1, p_2, \ldots, p_k} denote the individual (one- or two-sided) \eqn{p}-values of the \eqn{k} hypothesis tests to be combined, the combined \eqn{p}-value is then computed with \deqn{p_c = \min(1, \min(p_1, p_2, \ldots, p_k) \times k)}.

   The Bonferroni method does not assume that the \eqn{p}-values to be combined are independent. However, if the \eqn{p}-values are not independent, the method can become quite conservative (not reject often enough), depending on the correlation structure among the tests. In this case, one can adjust the method to account for such dependence (to bring the Type I error rate closer to some desired nominal significance level).

   \bold{Adjustment Based on the Effective Number of Tests}

   When \code{adjust} is set to \code{"nyholt"}, \code{"liji"}, \code{"gao"} or \code{"galwey"}, the Bonferroni method is adjusted based on an estimate of the effective number of tests (see \code{\link{meff}} for details on these methods for estimating the effective number of tests). In this case, argument \code{R} needs to be set to a correlation matrix that reflects the dependence structure among the tests.

   There is no general solution for constructing such a correlation matrix, as this depends on the type of test that generated the \eqn{p}-values and the sidedness of these tests. If the \eqn{p}-values are obtained from tests whose test statistics can be assumed to follow a multivariate normal distribution and a matrix is available that reflects the correlations among the test statistics, then the \code{\link{mvnconv}} function can be used to convert this correlation matrix into the correlations among the (one- or two-sided) \eqn{p}-values, which can then be passed to the \code{R} argument. See \sQuote{Examples}.

   Once the effective number of tests, \eqn{m}, is estimated based on \code{R} using one of the four methods described above, the combined \eqn{p}-value is then computed with \deqn{p_c = \min(1, \min(p_1, p_2, \ldots, p_k) \times m)}.

   Alternatively, one can also directly specify the effective number of tests via the \code{m} argument (e.g., if some other method not implemented in the \pkg{poolr} package is used to estimate the effective number of tests). Argument \code{R} is then irrelevant and doesn't need to be specified.

   \bold{Adjustment Based on an Empirically-Derived Null Distribution}

   When \code{adjust = "empirical"}, the combined \eqn{p}-value is computed based on an empirically-derived null distribution using pseudo replicates. This is appropriate if the test statistics that generated the \eqn{p}-values to be combined can be assumed to follow a multivariate normal distribution and a matrix is available that reflects the correlations among the test statistics (which is specified via the \code{R} argument). In this case, test statistics are repeatedly simulated from a multivariate normal distribution under the joint null hypothesis, converted into one- or two-sided \eqn{p}-values (depending on the \code{side} argument), and the Bonferroni method is applied. Repeating this process \code{size} times yields a null distribution based on which the combined \eqn{p}-value can be computed, or more precisely, estimated, since repeated applications of this method will yield (slightly) different results. To obtain a stable estimate of the combined \eqn{p}-value, \code{size} should be set to a large value (the default is \code{10000}, but this can be increased for a more precise estimate). If we consider the combined \eqn{p}-value an estimate of the \sQuote{true} combined \eqn{p}-value that would be obtained for a null distribution of infinite size, we can also construct a 95\% (pseudo) confidence interval based on a binomial distribution.

   If \code{batchsize} is unspecified, the null distribution is simulated in a single batch, which requires temporarily storing a matrix with dimensions \code{[size,k]}. When \code{size} and/or \code{k} is large, allocating the memory for this matrix might not be possible. Instead, one can specify a \code{batchsize} value, in which case a matrix with dimensions \code{[batchsize,k]} is repeatedly simulated until the desired size of the null distribution has been obtained.

   One can also specify a vector for the \code{size} argument, in which case one must also specify a corresponding vector for the \code{threshold} argument. In that case, a stepwise algorithm is used that proceeds as follows. For \code{j = 1, ..., length(size)},
   \enumerate{
      \item estimate the combined \eqn{p}-value based on \code{size[j]}
      \item if the combined \eqn{p}-value is \eqn{\ge} than \code{threshold[j]}, stop (and report the combined \eqn{p}-value), otherwise go back to 1.
   }
   By setting \code{size} to increasing values (e.g., \code{size = c(1000, 10000, 100000)}) and \code{threshold} to decreasing values (e.g., \code{threshold = c(.10, .01, 0)}), one can quickly obtain a fairly accurate estimate of the combined \eqn{p}-value if it is far from significant (e.g., \eqn{\ge} .10), but hone in on a more accurate estimate for a combined \eqn{p}-value that is closer to 0. Note that the last value of \code{threshold} should be 0 (and is forced to be inside of the function), so that the algorithm is guaranteed to terminate (hence, one can also leave out the last value of \code{threshold}, so \code{threshold = c(.10, .01)} would also work in the example above). One can also specify a single \code{threshold} (which is replicated as often as necessary depending on the length of \code{size}).

}
\value{
   An object of class \code{"poolr"}. The object is a list containing the following components:
   \item{p}{combined \eqn{p}-value.}
   \item{ci}{confidence interval for the combined \eqn{p}-value (only when \code{adjust = "empirical"}; otherwise \code{NULL}).}
   \item{k}{number of \eqn{p}-values that were combined.}
   \item{m}{estimate of the effective number of tests (only when \code{adjust} is one of \code{"nyholt"}, \code{"liji"}, \code{"gao"} or \code{"galwey"}; otherwise \code{NULL}).}
   \item{adjust}{chosen adjustment method.}
   \item{statistic}{value of the (adjusted) test statistic.}
   \item{fun}{name of calling function.}
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Bland, J. M., & Altman, D. G. (1995). Multiple significance tests: The Bonferroni method. \emph{British Medical Journal, 310}(6973), 170.
}
\examples{
# copy p-values and LD correlation matrix into p and r
# (see help(grid2ip) for details on these data)
p <- grid2ip.p
r <- grid2ip.ld

# apply the Bonferroni method
bonferroni(p)

# use mvnconv() to convert the LD correlation matrix into a matrix with the
# correlations among the (two-sided) p-values assuming that the test
# statistics follow a multivariate normal distribution with correlation
# matrix r (note: 'side = 2' by default in mvnconv())
mvnconv(r, target = "p", cov2cor = TRUE)[1:5,1:5] # show only rows/columns 1-5

# adjustment based on estimates of the effective number of tests
bonferroni(p, adjust = "nyholt", R = mvnconv(r, target = "p", cov2cor = TRUE))
bonferroni(p, adjust = "liji",   R = mvnconv(r, target = "p", cov2cor = TRUE))
bonferroni(p, adjust = "gao",    R = mvnconv(r, target = "p", cov2cor = TRUE))
bonferroni(p, adjust = "galwey", R = mvnconv(r, target = "p", cov2cor = TRUE))

# setting argument 'm' manually
bonferroni(p, m = 12)

# adjustment based on an empirically-derived null distribution (setting the
# seed for reproducibility)
set.seed(1234)
bonferroni(p, adjust = "empirical", R = r)

# generate the empirical distribution in batches of size 100
bonferroni(p, adjust = "empirical", R = r, batchsize = 100)

# using the stepwise algorithm
bonferroni(p, adjust = "empirical", R = r, size = c(1000, 10000, 100000),
           threshold = c(.10, .01))
}
\keyword{htest}
