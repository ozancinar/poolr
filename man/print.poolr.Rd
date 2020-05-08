\name{print.poolr}
\alias{print.poolr}
\title{Print Method for 'poolr' Objects}
\description{
   Print method for objects of class \code{"poolr"}.
}
\usage{
\method{print}{poolr}(x, digits=3, \dots)
}
\arguments{
   \item{x}{an object of class \code{"poolr"}.}
   \item{digits}{integer specifying the number of (significant) digits for rounding/presenting the results.}
   \item{\dots}{other arguments.}
}
\details{
   The output shows the combined \mjseqn{p}-value (with the specified number of significant digits), the test statistic (and its assumed null distribution), and the adjustment method that was applied.
}
\value{
   The function does not return an object.
}
\author{
   Ozan Cinar \email{ozancinar86@gmail.com} \cr
   Wolfgang Viechtbauer \email{wvb@wvbauer.com} \cr
}
\keyword{print}
