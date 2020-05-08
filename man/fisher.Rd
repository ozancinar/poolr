\name{fisher}
\alias{fisher}
\title{Fisher's Method}
\description{Function to carry out Fisher's method.\loadmathjax}
\usage{
fisher(p, adjust = "none", R, m, size = 10000, threshold, side = 2, batchsize, \dots)
}
\arguments{
   \item{p}{vector of length \mjseqn{k} with the (one- or two-sided) p-values to be combined.}
   \item{adjust}{character string to specify an adjustment method to account for dependence. The default is \code{"none"}, in which case no adjustment is applied. Methods \code{"nyholt"}, \code{"liji"}, \code{"gao"}, or \code{"galwey"} are adjustments based on an estimate of the effective number of tests (see \code{\link{meff}}). Adjustment method \code{"empirical"} uses an empirically-derived null distribution using pseudo replicates. Finally, method \code{"generalized"} uses a generalization of Fisher's method based on multivariate theory. See \sQuote{Details}.}
   \item{R}{a \mjeqn{k \times k}{k * k} symmetric matrix that reflects the dependence structure among the tests. Must be specified if \code{adjust} is set to something other than \code{"none"}. See \sQuote{Details}.}
   \item{m}{optional scalar (between 1 and \mjseqn{k}) to manually specify the effective number of tests (instead of estimating it via one of the methods described above).}
   \item{size}{size of the empirically-derived null distribution. Can also be a numeric vector of sizes, in which case a stepwise algorithm is used. This (and the following arguments) are only relevant when \code{adjust = "empirical"}.}
   \item{threshold}{numeric vector to specify the significance thresholds for the stepwise algorithm (only relevant when \code{size} is a vector).}
   \item{side}{scalar to specify the sidedness of the \mjseqn{p}-values that are used to simulate the null distribution (2, by default, for two-sided tests; 1 for one-sided tests).}
   \item{batchsize}{optional scalar to specify the batch size for generating the null distribution. When unspecified (the default), this is done in a single batch.}
   \item{\dots}{other arguments.}
}
\details{

   \bold{Fisher's Method}

   By default (i.e., when \code{adjust = "none"}), the function applies Fisher's method to the \mjseqn{p}-values (Fisher, 1932). Letting \mjseqn{p_1, p_2, \ldots, p_k} denote the individual (one- or two-sided) \mjseqn{p}-values of the \mjseqn{k} hypothesis tests to be combined, the test statistic is then computed with \mjdeqn{X^2 = -2 \sum_{i = 1}^k \ln(p_i).}{X^2 = -2 sum_{i=1}^k ln(p_i).} Under the joint null hypothesis, the test statistic follows a chi-square distribution with \mjseqn{2k} degrees of freedom which is used to compute the combined \mjseqn{p}-value.

   Fisher's method assumes that the \mjseqn{p}-values to be combined are independent. If this is not the case, the method can either be conservative (not reject often enough) or liberal (reject too often), depending on the correlation structure among the tests. In this case, one can adjust the method to account for such dependence (to bring the Type I error rate closer to some desired nominal significance level).

   \bold{Adjustment Based on the Effective Number of Tests}

   When \code{adjust} is set to \code{"nyholt"}, \code{"liji"}, \code{"gao"} or \code{"galwey"}, Fisher's method is adjusted based on an estimate of the effective number of tests (see \code{\link{meff}} for details on these methods for estimating the effective number of tests). In this case, argument \code{R} needs to be set to a correlation matrix that reflects the dependence structure among the tests.

   There is no general solution for constructing such a correlation matrix, as this depends on the type of test that generated the \mjseqn{p}-values and the sidedness of these tests. If the \mjseqn{p}-values are obtained from tests whose test statistics can be assumed to follow a multivariate normal distribution and a matrix is available that reflects the correlations among the test statistics, then the \code{\link{mvnconv}} function can be used to convert this correlation matrix into the correlations among the (one- or two-sided) \mjseqn{p}-values, which can then be passed to the \code{R} argument. See \sQuote{Examples}.

   Once the effective number of tests, \mjseqn{m}, is estimated based on \code{R} using one of the four methods described above, the test statistic of Fisher's method can be modified with \mjdeqn{\tilde{X}^2 = \frac{m}{k} \times X^2}{X'^2 = m/k * X^2} which is then assumed to follow a chi-square distribution with \mjseqn{2m} degrees of freedom.

   Alternatively, one can also directly specify the effective number of tests via the \code{m} argument (e.g., if some other method not implemented in the \pkg{poolr} package is used to estimate the effective number of tests). Argument \code{R} is then irrelevant and doesn't need to be specified.

   \bold{Adjustment Based on an Empirically-Derived Null Distribution}

   When \code{adjust = "empirical"}, the combined \mjseqn{p}-value is computed based on an empirically-derived null distribution using pseudo replicates. This is appropriate if the test statistics that generated the \mjseqn{p}-values to be combined can be assumed to follow a multivariate normal distribution and a matrix is available that reflects the correlations among the test statistics (which is specified via the \code{R} argument). In this case, test statistics are repeatedly simulated from a multivariate normal distribution under the joint null hypothesis, converted into one- or two-sided \mjseqn{p}-values (depending on the \code{side} argument), and Fisher's method is applied. Repeating this process \code{size} times yields a null distribution based on which the combined \mjseqn{p}-value can be computed, or more precisely, estimated, since repeated applications of this method will yield (slightly) different results. To obtain a stable estimate of the combined \mjseqn{p}-value, \code{size} should be set to a large value (the default is \code{10000}, but this can be increased for a more precise estimate). If we consider the combined \mjseqn{p}-value an estimate of the \sQuote{true} combined \mjseqn{p}-value that would be obtained for a null distribution of infinite size, we can also construct a 95\% (pseudo) confidence interval based on a binomial distribution.

   If \code{batchsize} is unspecified, the null distribution is simulated in a single batch, which requires temporarily storing a matrix with dimensions \code{[size,k]}. When \code{size} and/or \code{k} is large, allocating the memory for this matrix might not be possible. Instead, one can specify a \code{batchsize} value, in which case a matrix with dimensions \code{[batchsize,k]} is repeatedly simulated until the desired size of the null distribution has been obtained.

   One can also specify a vector for the \code{size} argument, in which case one must also specify a corresponding vector for the \code{threshold} argument. In that case, a stepwise algorithm is used that proceeds as follows. For \code{j = 1, ..., length(size)},
   \enumerate{
      \item estimate the combined \mjseqn{p}-value based on \code{size[j]}
      \item if the combined \mjseqn{p}-value is \mjseqn{\ge} than \code{threshold[j]}, stop (and report the combined \mjseqn{p}-value), otherwise go back to 1.
   }
   By setting \code{size} to increasing values (e.g., \code{size = c(1000, 10000, 100000)}) and \code{threshold} to decreasing values (e.g., \code{threshold = c(.10, .01, 0)}), one can quickly obtain a fairly accurate estimate of the combined \mjseqn{p}-value if it is far from significant (e.g., \mjseqn{\ge} .10), but hone in on a more accurate estimate for a combined \mjseqn{p}-value that is closer to 0. Note that the last value of \code{threshold} should be 0 (and is forced to be inside of the function), so that the algorithm is guaranteed to terminate (hence, one can also leave out the last value of \code{threshold}, so \code{threshold = c(.10, .01)} would also work in the example above). One can also specify a single \code{threshold} (which is replicated as often as necessary depending on the length of \code{size}).

   \bold{Adjustment Based on Multivariate Theory}

   When \code{adjust = "generalized"}, Fisher's method is computed based on a Satterthwaite approximation that accounts for the dependence among the tests, assuming that the test statistics that generated the \mjseqn{p}-values follow a multivariate normal distribution. In that case, \code{R} needs to be set equal to a matrix that contains the covariances among the \mjeqn{-2 \ln(p_i)}{-2 ln(p_i)} values. If a matrix is available that reflects the correlations among the test statistics, this can be converted into the required covariance matrix using the \code{\link{mvnconv}} function. See \sQuote{Examples}.

   This generalization of Fisher's method is sometimes called Brown's method, based on Brown (1975), although the paper only describes the method for combining one-sided p-values. Both one- and two-sided versions of Brown's method are implemented in \pkg{poolr}.

}
\value{
   An object of class \code{"poolr"}. The object is a list containing the following components:
   \item{p}{combined \mjseqn{p}-value.}
   \item{ci}{confidence interval for the combined \mjseqn{p}-value (only when \code{adjust = "empirical"}; otherwise \code{NULL}).}
   \item{k}{number of \mjseqn{p}-values that were combined.}
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
   Brown, M. B. (1975). 400: A method for combining non-independent, one-sided tests of significance. \emph{Biometrics, 31}(4), 987-992.

   Fisher, R. A. (1932). \emph{Statistical Methods for Research Workers} (4th ed.). Edinburgh: Oliver and Boyd.
}
\examples{
# copy p-values and LD correlation matrix into p and r
# (see help(grid2ip) for details on these data)
p <- grid2ip.p
r <- grid2ip.ld

# apply Fisher's method
fisher(p)

# use mvnconv() to convert the LD correlation matrix into a matrix with the
# correlations among the (two-sided) p-values assuming that the test
# statistics follow a multivariate normal distribution with correlation
# matrix r (note: 'side = 2' by default in mvnconv())
mvnconv(r, target = "p", cov2cor = TRUE)[1:5,1:5] # show only rows/columns 1-5

# adjustment based on estimates of the effective number of tests
fisher(p, adjust = "nyholt", R = mvnconv(r, target = "p", cov2cor = TRUE))
fisher(p, adjust = "liji",   R = mvnconv(r, target = "p", cov2cor = TRUE))
fisher(p, adjust = "gao",    R = mvnconv(r, target = "p", cov2cor = TRUE))
fisher(p, adjust = "galwey", R = mvnconv(r, target = "p", cov2cor = TRUE))

# setting argument 'm' manually
fisher(p, m = 12)

# adjustment based on an empirically-derived null distribution (setting the
# seed for reproducibility)
set.seed(1234)
fisher(p, adjust = "empirical", R = r)

# generate the empirical distribution in batches of size 100
fisher(p, adjust = "empirical", R = r, batchsize = 100)

# using the stepwise algorithm
fisher(p, adjust = "empirical", R = r, size = c(1000, 10000, 100000), threshold = c(.10, .01))

# use mvnconv() to convert the LD correlation matrix into a matrix with the
# covariances among the (two-sided) '-2ln(p_i)' values assuming that the
# test statistics follow a multivariate normal distribution with correlation
# matrix r (note: 'side = 2' by default in mvnconv())
mvnconv(r, target = "m2lp")[1:5,1:5] # show only rows/columns 1-5

# adjustment based on generalized method
fisher(p, adjust = "generalized", R = mvnconv(r, target = "m2lp"))

# when using mvnconv() inside fisher() with adjust = "generalized", the
# 'target' argument is automatically set and doesn't need to be specified
fisher(p, adjust = "generalized", R = mvnconv(r))
}
\keyword{htest}
