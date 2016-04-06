\name{one.pairscan.parallel}
\alias{one.pairscan.parallel}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Perform regression analysis for one phenotype and all pairs of markers using a kinship correction.
}
\description{
	This is an internal function that performs a single pairwise scan. It is used for individual pair scans as well as permutation testing.
}
\usage{
one.pairscan.parallel(data.obj, phenotype.vector, genotype.matrix, 
int = NULL, covar.vector = NULL, paired.markers, 
n.perm = 0, run.parallel = TRUE, verbose = FALSE, n.cores = 2)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
  \item{phenotype.vector}{
	A vector of phenotype values, one entry for each individual.
}
  \item{genotype.matrix}{
	A matrix of genotype values with individuals in rows and markers in columns. Matrix entries contain the probability of the reference allele at each position for each individual.
}
\item{int}{
	The interaction term used in the kinship correction.
	}
  \item{covar.vector}{
	A vector of 1s and 0s indicating which markers should be used as covariates (1) and which should not (0),
	}
	\item{paired.markers}{
		A numeric matrix with two columns specifying which marker pairs should be scanned.
	}
  \item{n.perm}{
	The number of permutations to perform. May be 0.	
	}
	\item{run.parallel}{
	A logical value indicating whether to use parallel processing.	
	}
  \item{verbose}{
	A logical value indicating whether the progress of the scan should be printed to the screen.
}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing.	}
}
\value{
This function returns a list with two elements:
\item{pairscan.results}{The results of the pairwise scan on the provided phenotype and genotypes.}
\item{pairscan.perm}{The results of the permutations of the pairwise scan on the provided phenotype and genotypes.}

Each of the results elements itself a list of 3 elements:
\item{pairscan.effects}{A table of effects of each marker pair. The columns list the effects in the following order: marker1, marker2, the variance of marker1, the covariance of marker1 and marker2, the variance of marker2, the covariance of marker1 and the interaction effect, the covariance between marker2 and ther interaction effect, and the variance of the interaction.}
\item{pairscan.se}{A table of the standard errors from the test on each marker pair. The columns are identical to those described for pairscan.effects}
\item{model.covariance}{This is a table in which each row is the linearized matrix of the variance-covariance matrix of each pairwise regression.}

The results element for the permutation tests has the same structure as for the pairwise scan except that each row represents the results of one permutation.

}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010

Kang, H.M., Zaitlen, N.A., Wade, C.M., Kirby, A., Heckerman, D., Daly, M. J., and Eskin, E. (2008). Efficitne control of population structure in model organism association mapping. Genetics, 178(3), 1709-1723.
}
\seealso{
\code{\link{pairscan}}
}

% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}
\keyword{regression}% __ONLY ONE__ keyword per line
