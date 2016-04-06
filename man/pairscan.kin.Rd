\name{pairscan.kin}
\alias{pairscan.kin}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Performs a pairscan with kinship correction.
	}
\description{
This function performs the pairwise regression on all selected marker pairs. The phenotypes used can be either eigentraits or raw phenotypes.
}
\usage{
pairscan.kin(data.obj, geno.obj = NULL, covar, 
scan.what, marker.pairs, kin.full.geno, 
sample.kinship, num.kin.samples, n.per.sample, 
verbose = TRUE, n.cores = 2)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{geno.obj}{
	The object in which the genotype matrix and marker information are stored. See \code{\link{read.geno}}.
	}
	\item{covar}{
		A vector of character strings indicating which covariates to use.
	}
  \item{scan.what}{
	A character string uniquely identifying whether eigentraits or raw traits should be scanned.
	}
	\item{marker.pairs}{
		A two-column matrix indicating which marker pairs to test in the pairscan.
	}
	\item{kin.full.geno}{
		A logical value indicating whether a kinship matrix for the entire genome should be calculated. This is TRUE if phenotypic covariates are being tested.
	}
	\item{sample.kinship}{
		A logical value indicating whether kinship matrices should be sampled (TRUE) or calculated directly (FALSE). Sampling is faster when very large genotype matrices are being used.
	}
	\item{num.kin.samples}{
		An integer indicating how many samples to use in calculating the kinship matrices if sample.kinship is TRUE.
	}
	\item{n.per.sample}{
		An integer indicating how many markers should be sampled in calculating the kinship matrices if sample.kinship is TRUE.
	}
  \item{verbose}{
	A logical value indicating whether the progress of the scan should be printed to the screen.
	}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing.	}
}
\value{
This function adds an element to the pairscan object reporting the results of the pair-wise scan:
\item{pairscan.effects}{A matrix with five columns indicating the names of markers 1 and 2, their effect sizes, and the effect size of their interaction}
\item{pairscan.se}{A table of the standard errors from the test on each marker pair. The columns are identical to those described for pairscan.effects}
\item{model.covariance}{This is a table in which each row is the linearized matrix of the variance-covariance matrix of each pairwise regression.}
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{pairscan}}, \code{\link{pairscan.noKin}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

