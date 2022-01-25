\name{poolr-package}
\alias{poolr-package}
\alias{poolr}
\docType{package}
\title{Methods for Pooling P-Values from (Dependent) Tests}
\description{
   \loadmathjax The \pkg{poolr} package contains functions for pooling/combining the results (i.e., \mjseqn{p}-values) from (dependent) hypothesis tests. Included are Fisher's method, Stouffer's method, the inverse chi-square method, the Bonferroni method, Tippett's method, and the binomial test. Each method can be adjusted based on an estimate of the effective number of tests or using empirically-derived null distribution using pseudo replicates. For Fisher's, Stouffer's, and the inverse chi-square method, direct generalizations based on multivariate theory are also available (leading to Brown's method, Strube's method, and the generalized inverse chi-square method). For more details, see:

   \itemize{
   \item \code{\link{fisher}}: for Fisher's method (and Brown's method)
   \item \code{\link{stouffer}}: for Stouffer's method (and Strube's method)
   \item \code{\link{invchisq}}: for the inverse chi-square method
   \item \code{\link{bonferroni}}: for the Bonferroni method
   \item \code{\link{tippett}}: for Tippett's method
   \item \code{\link{binomtest}}: for the binomial test
   }

   Note that you can also read the documentation of the package online at \url{https://ozancinar.github.io/poolr/} (where it is nicely formatted and the output from all examples is provided).
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Brown, M. B. (1975). 400: A method for combining non-independent, one-sided tests of significance. \emph{Biometrics, 31}(4), 987--992.

   Cinar, O. & Viechtbauer, W. (2022). The poolr package for combining independent and dependent p values. \emph{Journal of Statistical Software}, \bold{101}(1), 1--42. \verb{https://doi.org/10.18637/jss.v101.i01}

   Fisher, R. A. (1932). \emph{Statistical Methods for Research Workers} (4th ed.). Edinburgh: Oliver and Boyd.

   Lancaster, H. O. (1961). The combination of probabilities: An application of orthonormal functions. \emph{Australian Journal of Statistics, 3}(1), 20--33.

   Strube, M. J. (1985). Combining and comparing significance levels from nonindependent hypothesis tests. \emph{Psychological Bulletin, 97}(2), 334--341.

   Tippett, L. H. C. (1931). \emph{Methods of Statistics}. London: Williams Norgate.

   Wilkinson, B. (1951). A statistical consideration in psychological research. \emph{Psychological Bulletin, 48}(2), 156--158.
}
\keyword{package}
