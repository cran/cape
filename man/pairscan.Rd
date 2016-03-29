\name{pairscan}
\alias{pairscan}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Perform regressions for all pairs of markers and all phenotypes.
}
\description{
This function performs the pairwise regression on all selected marker pairs. The phenotypes used can be either eigentraits or raw phenotypes. Permutation testing is also performed, and kinship corrections are implemented if requested.
}
\usage{
pairscan(data.obj, geno.obj = NULL, covar = NULL, 
scan.what = c("eigentraits", "raw.traits"), total.perm = NULL, 
min.per.genotype = NULL, max.pair.cor = NULL, 
n.top.markers = NULL, use.kinship = FALSE, kin.full.geno = TRUE, 
sample.kinship = TRUE, num.kin.samples = 100, n.per.sample = 100, 
verbose = FALSE, num.pairs.limit = 1e6, num.perm.limit = 1e7,
overwrite.alert = TRUE, n.cores = NULL)
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
		A vector of character strings indicating which covariates should be used in the regressions.
	}
  \item{scan.what}{
	A character string uniquely identifying whether eigentraits or raw traits should be scanned.
}
  \item{total.perm}{
	The total number of permutations to be performed.
}
  \item{min.per.genotype}{
	The minimum number of individuals allowable per genotype. If for a given marker pair, one of the genotypes is underrepresented, the marker pair is not tested. If this value is NULL, max.pair.cor must have a numeric value.
}
	\item{max.pair.cor}{
		A numeric value between 0 and 1 indicating the maximum Pearson correlation that two markers are allowed. If the correlation between a pair of markers exceeds this threshold, the pair is not tested. If this value is set to NULL, min.per.genotype must have a numeric value.
	}
	\item{n.top.markers}{
	This optional integer is used in the generation of a null distribution. To generate the null, a singlescan is performed on permuted data and the top ranking markers are used in a pairscan. n.top.markers specifies how many of these markers should be chosen. If NULL, n.top.markers is the number of markers being scanned in the pairscan.
	}
	\item{use.kinship}{
		A logical value indicating whether a kinship correction should be implemented.
		}
	\item{kin.full.geno}{
		A logical value indicating whether the entire genotype matrix should be used to calculate kinship corrections. If FALSE, only the subset of genetic markers being scanned will be used to calculate kinship corrections.
	}
	\item{sample.kinship}{
		A logical value indicating whether kinship matrices should be sampled (TRUE) or calculated directly (FALSE). Sampling kinshp matrices can be faster in the case of very large genotype matrices.
	}
	\item{num.kin.samples}{
		If sample.kinship is TRUE, this integer indicates how many samples should be used to calculate kinship matrices.
	}
	\item{n.per.sample}{
		If sample.kinship is TRUE, this integer indicates how many markers should be used in each sample of the genotype matrix for calculating kinship matrices.
	}
  \item{verbose}{
	A logical value indicating whether the progress of the scan should be printed to the screen.
	}
\item{num.pairs.limit}{
	A number indicating the maximum number of pairs to scan. If the number of pairs exceeds this threshold, the function asks for confirmation before proceeding with the pairwise scan.
	}
\item{num.perm.limit}{
	A number indicating the maximum number of total permutations that will be performed. If the number of total permutations exceeds this threshold, the function asks for confirmation before proceeding with the pairwise scan.
	}
	\item{overwrite.alert}{
	If TRUE, this triggers an alert warning the user that the output of this function should be saved separately from the data.obj.
	}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing. If NULL, the choice is made automatically.
	}
}
\details{
Not all marker pairs are necessarily tested. Before markers are tested for interaction, they are checked for several conditions. Pairs are discarded if (1) at least one of the markers is on the X chromosome, or (2) there are fewer than min.per.genotype individuals in any of the genotype combinations.
}
\value{
IMPORTANT NOTE: Prior versions of this function modified the data.obj. The new version creates a separate pairscan.obj, and the output of this function should NOT be assigned to the data.obj.
The pairscan.obj contains the following elements:

\item{pairscan.results}{The results of the pairwise scan on the provided phenotype and genotypes.}
If permutations have been performed (n.perm > 0), an additional element is added to the object reporting the results of the permutation tests:
\item{pairscan.perm}{The results of the permutations of the pairwise scan on the provided phenotype and genotypes.}
Each of these results elements is itself a list of 3 elements:
\item{pairscan.effects}{A table of effects of each marker pair. The columns list the effects in the following order: marker1, marker2, the variance of marker1, the covariance of marker1 and marker2, the variance of marker2, the covariance of marker1 and the interaction effect, the covariance between marker2 and ther interaction effect, and the variance of the interaction.}
\item{pairscan.se}{A table of the standard errors from the test on each marker pair. The columns are identical to those described for pairscan.effects}
\item{model.covariance}{This is a table in which each row is the linearized matrix of the variance-covariance matrix of each pairwise regression.}
\item{pairs.tested}{
	A two-column matrix identifying the marker pairs tested for each permutation.
	}

The results element for the permutation tests has the same structure as for the pairwise scan except that each row represents the results of one permutation.
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{select.markers.for.pairscan}}, \code{\link{plotPairscan}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{regression}

