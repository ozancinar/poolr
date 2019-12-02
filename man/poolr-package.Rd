\name{poolr-package}
\alias{poolr-package}
\alias{poolr}
\docType{package}
\title{Methods for Pooling P-Values from (Dependent) Tests}
\description{
   The \pkg{poolr} package is a collection of functions for pooling/combining the results (i.e., p-values) from (dependent) hypothesis tests. Included are Fisher's method, Stouffer's method, the inverse chi-squared method, the Bonferroni method, Tippett's method, and the binomial test. Each method can be adjusted based on an estimate of the effective number of tests or using empirically derived null distribution using pseudo replicates. For Fisher's, Stouffer's, and the inverse chi-squared method, direct generalizations based on multivariate theory are also available (leading to Brown's method, Strube's method, and the generalized inverse chi-squared method).
}
\section{Purpose}{
   The problem to be addressed can be succinctly stated as follows. Let \eqn{p_i} denote the \eqn{p}-value corresponding to the \eqn{i}th null hypothesis, \eqn{H_i}, where \eqn{i = 1, \ldots, k}. The tests being used to generate the \eqn{p_i} values are assumed to have nominal properties (i.e., the Type I error rates are equal to the chosen \eqn{\alpha} value), which implies \eqn{p_i \sim Uniform(0,1)} when \eqn{H_i} is true. The goal is to combine the \eqn{p}-values into a single \eqn{p}-value to test the joint null hypothesis that \eqn{H_i} is true for all \eqn{k} tests.
}
\section{Methods}{
   Blah.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Bland, J. M., & Altman, D. G. (1995). \emph{Multiple significance tests: The Bonferroni method}. British Medical Journal, 310(6973), 170.

   Brown, MB (1975): \emph{400: A method for combining non-independent, one-sided tests of significance}. Biometrics, 987-992.

   Fisher, R. A. (1932). \emph{Statistical Methods for Research Workers} (4th ed.). Edinburgh: Oliver and Boyd.

   Mehta, C. R., & Patel, N. R. (1998). \emph{Exact inference for categorical data}. Encyclopedia of Biostatistics, 2, 1411-1422.

   Strube, M. J. (1985). \emph{Combining and comparing significance levels from nonindependent hypothesis tests}. Psychological Bulletin, 97(2), 334.

   Tippett, L. H. C. (1931). \emph{Methods of Statistics}. London: Williams Norgate.

   Wilkinson, B. (1951). \emph{A statistical consideration in psychological research}. Psychological Bulletin, 48(2), 156-158.
}
\keyword{package}
