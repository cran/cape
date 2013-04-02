\name{pairscan}
\alias{pairscan}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Perform regressions for all pairs of markers and all phenotypes.
}
\description{
This function performs the pairwise regression on all selected marker pairs. The phenotypes used can be either eigentraits or raw phenotypes. Permutation testing is also performed.
}
\usage{
pairscan(data.obj, scan.what = c("eigentraits", "raw.traits"), 
n.perm = NULL, min.per.genotype = 6, verbose = FALSE, 
num.pairs.limit = 1e4, num.perm.limit = 1e7)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
  \item{scan.what}{
	A character string uniquely identifying whether eigentraits or raw traits should be scanned.
}
  \item{n.perm}{
	The number of permutations.
}
  \item{min.per.genotype}{
	The minimum number of individuals allowable per genotype. If for a given marker pair, one of the genotypes is underrepresented, the marker pair is not tested.
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
}
\details{
Not all marker pairs are necessarily tested. Before markers are tested for interaction, they are checked for several conditions. Pairs are discarded if (1) at least one of the markers is on the X chromosome, or (2) there are fewer than min.per.genotype individuals in any of the genotype combinations.
}
\value{
This function adds an element to the data object reporting the results of the pair-wise scan:
\item{pairscan.results}{The results of the pairwise scan on the provided phenotype and genotypes.}

If permutations have been performed (n.perm > 0), an additional element is added to the object reporting the results of the permutation tests:
\item{pairscan.perm}{The results of the permutations of the pairwise scan on the provided phenotype and genotypes.}

Each of these results elements is itself a list of 3 elements:
\item{pairscan.effects}{A table of effects of each marker pair. The columns list the effects in the following order: marker1, marker2, the variance of marker1, the covariance of marker1 and marker2, the variance of marker2, the covariance of marker1 and the interaction effect, the covariance between marker2 and ther interaction effect, and the variance of the interaction.}
\item{pairscan.se}{A table of the standard errors from the test on each marker pair. The columns are identical to those described for pairscan.effects}
\item{model.covariance}{This is a table in which each row is the linearized matrix of the variance-covariance matrix of each pairwise regression.}


The results element for the permutation tests has the same structure as for the pairwise scan except that each row represents the results of one permutation.
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{select.markers.for.pairscan}}, \code{\link{plotPairscan}}
}
\examples{
# not run
# obesity.cross <- pairscan(obesity.cross, scan.what = "eigentraits", n.perm = 10, 
# min.per.genotype = 6, use.pairs.threshold = TRUE, verbose = TRUE)
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{regression}

