\name{pairscan.null}
\alias{pairscan.null}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Generate a null distribution for the pairscan.
}
\description{
This function performs the pairwise regression on all permuted pairs of markers to generate a null distribution. The marker pairs are selected by permutation. A singlescan is performed with a permuted phenotype, and the N markers with the highest effect sizes are selected. N is equal to the number of markers being tested in the pairscan. The selected markers are then tested in pairs for association with the permuted phenotype, and the results are collected to form a null distribution.
}
\usage{
pairscan.null(data.obj, geno.mat, covar = NULL, 
scan.what = c("eigentraits", "raw.traits"), 
total.perm = NULL, n.top.markers = 50, 
max.pair.cor, min.per.genotype, verbose = FALSE,
n.cores = NULL)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{geno.mat}{
	A matrix of genotype values with individuals in rows and markers in columns. Matrix entries contain the probability of the reference allele at each position for each individual.
	}
  \item{scan.what}{
	A character string uniquely identifying whether eigentraits or raw traits should be scanned.
}
  \item{total.perm}{
	The total number of permutations to be run across all marker pairs.
}
	\item{top.n.markers}{
		The number of markers to be selected for permutation testing. This is the same number as selected for the pairscan.
	}
	 \item{min.per.genotype}{
	The minimum number of individuals allowable per genotype. If for a given marker pair, one of the genotypes is underrepresented, the marker pair is not tested. If this value is NULL, max.pair.cor must have a numeric value.
}
	\item{max.pair.cor}{
		A numeric value between 0 and 1 indicating the maximum Pearson correlation that two markers are allowed. If the correlation between a pair of markers exceeds this threshold, the pair is not tested. If this value is set to NULL, min.per.genotype must have a numeric value.
	}
  \item{verbose}{
	A logical value indicating whether the progress of the scan should be printed to the screen.
}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing. If NULL, the choice is made automatically.
	}
}
\value{
This function adds an element to the pairscan object reporting the results of the pair-wise permutations:

\item{pairscan.perm}{The results of the permutations of the pairwise scan on the provided phenotype and genotypes.}

Each of these results elements is itself a list of 3 elements:
\item{pairscan.effects}{A table of effects of each marker pair. The columns list the effects in the following order: marker1, marker2, the variance of marker1, the covariance of marker1 and marker2, the variance of marker2, the covariance of marker1 and the interaction effect, the covariance between marker2 and ther interaction effect, and the variance of the interaction.}
\item{pairscan.se}{A table of the standard errors from the test on each marker pair. The columns are identical to those described for pairscan.effects}
\item{model.covariance}{This is a table in which each row is the linearized matrix of the variance-covariance matrix of each pairwise regression.}

The results element for the permutation tests has the same structure as for the pairwise scan except that each row represents the results of one permutation.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{pairscan}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

