\name{mvnlookup}
\docType{data}
\alias{mvnlookup}
\title{Lookup Tables for MVN Test Statistics}
\description{Lookup tables that include the covariances between two test statistics that follow a distribution assumed by the method of combining p-values respective for the correlation between the test statistics that follow multivariate normal distribution.}

\usage{
	data(mvnlookup)
}

\format{The data frame contains the following columns:
\tabular{lll}{
	\bold{rhos}     	\tab \code{numeric}  	\tab correlation among the mvn test statistics \cr
	\bold{m2lp_1}   	\tab \code{numeric} 	\tab covariances among -2log(p) (one-sided test statistics) \cr
	\bold{m2lp_2}   	\tab \code{numeric} 	\tab covariances among -2log(p) (two-sided test statistics) \cr
	\bold{z_1}   		\tab \code{numeric} 	\tab covariances among z (one-sided test statistics) \cr
	\bold{z_2}   		\tab \code{numeric} 	\tab covariances among z (two-sided test statistics) \cr
	\bold{chisq1_1}   	\tab \code{numeric} 	\tab covariances among Chi-Square (1) (one-sided test statistics) \cr
	\bold{chisq1_2}   	\tab \code{numeric} 	\tab covariances among Chi-Square (1) (two-sided test statistics) \cr
}
}
\details{
   Lookup tables for mvn test statistics.
}
\keyword{datasets}
