\name{fisher}
\alias{fisher}
\title{Fisher's Method}
\description{Function to carry out Fisher's method.}
\usage{
fisher(p, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to apply Fisher's method.
}
\value{
   The function returns the (pooled) p-value based on Fisher's method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Fisher, R. A. (1932). \emph{Statistical methods for research workers} (4th ed.). Oliver and Boyd.
}
\examples{
p <- c(.013, .022, .124, .278)
fisher(p)
}
\keyword{htest}
