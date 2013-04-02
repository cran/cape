\name{cape-package}
\alias{cape-package}
\alias{cape}
\docType{package}
\title{ 
	Combinatorial Analysis of Epistasis and Pleiotropy
}
\description{
This package infers predictive networks between genetic variants and between genetic variants and phenotypes. It uses complementary information of pleiotropic gene variants across different phenotypes to resolve models of epistatic interactions between genetic variants. To do this, cape reparameterizes main effect and interaction coefficients from a pairwise variant regressions into directed influence parameters. These parameters describe how gene variants influence each other, in terms of suppression and enhancement, as well as how gene variants influence phenotypes.
}
\details{
\tabular{ll}{
Package: \tab cape\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-12-30\cr
License: \tab GPL-3\cr
}
The cape analysis begins by reading in a genetic data set with \code{\link{read.population}}. The data are converted into a data object, to which results are added throughout the analysis. This data object is an argument in most functions, and is referred to as \code{data.obj}. Because this package uses pleiotropy to resolve models of epistasis, the phenotypes used should have common underlying molecular players, but not be perfectly correlated. In general phenotypes should be correlated with a Pearson r between 0.4 and 0.8. The phenotypes of interest are then decomposed into eigentraits using \code{\link{get.eigentraits}}. Any number of phenotypes can be decomposed, but \code{cape} requires between two and 12 eigentraits for the analysis. The phenotype decomposition into eigentraits maximizes the complementary information between the phenotypes. Before investigating epistatic interactions through a pair-wise scan of the genetic variants, a single-variant scan (\code{\link{singlescan}}) is run. This allows for thresholding of markers for the pair scan if the cross is prohibitively large to test all pairs of variants. The single-variant scan also allows selection of variants with very large main effects to be used as covariates in the pair scan. The pair scan (\code{\link{pairscan}}) performs a regression on each variant pair. Finally, the coefficients from the pairwise scan are reparameterized to yield directional influences between variants and from variants to phenotypes. The final result is an asymmetric adjacency matrix describing these variant influences. The p values of these influences are corrected for multiple tests.
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
Maintainer: Anna L. Tyler <anna.tyler@jax.org>
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010}

\keyword{ package }
\examples{
data(obesity.cross)
str(obesity.cross)
}
\keyword{package}
