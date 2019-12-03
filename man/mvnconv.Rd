\name{mvnconv}
\alias{mvnconv}
\title{Convert Correlations Among Multivariate Normal Test Statistics to Covariances for Various Target Statistics}
\description{Function to convert a matrix with the correlations among multivariate normal test statistics to a matrix with the covariances among various target statistics.}
\usage{
mvnconv(R, side = 2, target, cov2cor = FALSE)
}
\arguments{
   \item{R}{a \eqn{k \times k} symmetric matrix that contains the correlations among the test statistics.}
   \item{side}{scalar to specify the sidedness of the \eqn{p}-values that are obtained from the test statistics (2, by default, for two-sided tests; 1 for one-sided tests).}
   \item{target}{the target statistic for which the covariances are calculated (either \code{"p"}, \code{"m2lp"}, \code{"chisq1"}, or \code{"z"}). See \sQuote{Details}.}
   \item{cov2cor}{logical to indicate whether to convert the covariance matrix to a correlation matrix (default is \code{FALSE}).}
}
\details{

   The function converts a matrix with the correlations among multivariate normal test statistics to a matrix with the covariances among various target statistics. In particular, assume

   \if{html}{\figure{mvn.png}{options: width=250}}
   \if{latex}{\figure{mvn.png}{options: width=2in}}
   \if{text}{[t_i, t_j]' ~ MVN([0,0]', [1, rho_ij | rho_ij, 1])}

   is the joint distribution for test statistics \eqn{t_i} and \eqn{t_j}. For \code{side = 1}, let \eqn{p_i = 1 - \Phi(t_i)} and \eqn{p_j = 1 - \Phi(t_j)} where \eqn{\Phi(\cdot)} denotes the cumulative distribution function of a standard normal distribution. For \code{side = 2}, let \eqn{p_i = 2(1 - \Phi(|t_i|))} and \eqn{p_j = 2(1 - \Phi(|t_j|))}. These are simply the one- and two-sided \eqn{p}-values corresponding to \eqn{t_i} and \eqn{t_j}.

   If \code{target="p"}, the function computes \eqn{\mbox{Cov}[p_i, p_j]}.

   If \code{target="m2lp"}, the function computes \eqn{\mbox{Cov}[-2 \ln(p_i), -2 \ln(p_j)]}.

   If \code{target="chisq1"}, the function computes \eqn{\mbox{Cov}[F^{-1}(1 - p_i, 1), F^{-1}(1 - p_j, 1)]}, where \eqn{F^{-1}(\cdot,1)} denotes the inverse of the cumulative distribution function of a chi-square distribution with one degree of freedom.

   If \code{target="z"}, the function computes \eqn{\mbox{Cov}[\Phi^{-1}(1 - p_i), \Phi^{-1}(1 - p_j)]}, where \eqn{\Phi^{-1}(\cdot)} denotes the inverse of the cumulative distribution function of a standard normal distribution.

}
\value{
   The function returns the covariance matrix (or the correlation matrix if \code{cov2cor=TRUE}).
}
\note{
   Since computation of the covariances requires numerical integration, the function doesn't actually compute these covariances on the fly. Instead, it uses the \code{\link{mvnlookup}} lookup table, which contains the covariances.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\examples{
# illustrative correlation matrix
R <- matrix(c(   1,  0.8,  0.5, -0.3,
               0.8,    1,  0.2,  0.4,
               0.5,  0.2,    1, -0.7,
              -0.3,  0.4, -0.7,    1), nrow = 4, ncol = 4)

# convert R into covariance matrices for the chosen targets
mvnconv(R, target = "p")
mvnconv(R, target = "m2lp")
mvnconv(R, target = "chisq1")
mvnconv(R, target = "z")

# convert R into correlation matrices for the chosen targets
mvnconv(R, target = "p",      cov2cor = TRUE)
mvnconv(R, target = "m2lp",   cov2cor = TRUE)
mvnconv(R, target = "chisq1", cov2cor = TRUE)
mvnconv(R, target = "z",      cov2cor = TRUE)
}
\keyword{htest}
