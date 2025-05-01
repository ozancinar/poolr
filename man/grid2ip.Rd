\name{grid2ip}
\docType{data}
\alias{grid2ip.p}
\alias{grid2ip.ld}
\alias{grid2ip.geno}
\alias{grid2ip.pheno}
\alias{grid2ip}
\title{Results from testing the association between depressive symptoms and 23 SNPs in the GRID2IP gene}
\description{Results from testing the association between depressive symptoms (as measured with the CES-D scale) and 23 single-nucleotide polymorphisms (SNPs) in the GRID2IP gene based on a sample of 886 adolescents (Van Assche et al., 2017).}
\usage{grid2ip.p
grid2ip.ld
grid2ip.geno
grid2ip.pheno}
\format{Object \code{grid2ip.p} is a vector with the 23 \eqn{p}-values of the tests (two-sided). Object \code{grid2ip.ld} contains a matrix with the linkage disequilibrium (LD) correlations among the 23 SNPs. \code{grid2ip.geno} is a matrix that contains the genotypes of the adoloscents for the 23 SNPs. \code{grid2ip.pheno} is a vector with the phenotype for the adoloscents (the log-transformed CES-D scale values).}
\references{
   Van Assche, E., Moons, T., Cinar, O., Viechtbauer, W., Oldehinkel, A. J., Van Leeuwen, K., Verschueren, K., Colpin, H., Lambrechts, D., Van den Noortgate, W., Goossens, L., Claes, S., & van Winkel, R. (2017). Gene-based interaction analysis shows GABAergic genes interacting with parenting in adolescent depressive symptoms. \emph{Journal of Child Psychology and Psychiatry, 58}(12), 1301--1309. \verb{https://doi.org/10.1111/jcpp.12766}
}
\keyword{datasets}
