\name{minp}
\alias{minp}
\title{Minimum p}
\description{Function to combine the p-values with the strongest evidence, i.e. the smallest p-values.}
\usage{
minp(p, adjust, \dots)
}
\arguments{
   \item{p}{vector with the (one- or two-sided) p-values.}
   \item{adjust}{adjustment method while combining the p-values. No adjustments are done to consider the dependencies between the p-values, as default.}
   \item{\dots}{other arguments.}
}
\details{
   The function can be used to combine a set of p-values with the strongest evidence which is determined by the smallest p-value.
}
\value{
   The function returns the (pooled) p-value based on minimum p-value method.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\references{
   Tippett, L. H. C. (1931). \emph{The Methods of Statistics}. The methods of statistics.
}
\examples{
p <- c(.013, .022, .124, .278)
minp(p = p)
}
\keyword{htest}
