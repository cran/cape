\name{select.markers.for.pairscan}
\alias{select.markers.for.pairscan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
A required step that filters variable and non-redundant markers for the pairscan
}
\description{
This function selects markers for the pair scan. If the markers are to be thresholded by a significance cutoff, this function 
filters them. It then checks to make sure all pairs of markers have had at least one recombination between them, and that all 
markers are variable across individuals. Mono-allelic markers are removed. If any pair of markers carry identical genotype 
information, the first marker of the pair is discarded.
}
\usage{
select.markers.for.pairscan(data.obj, use.pairs.threshold = FALSE, 
specific.markers = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{use.pairs.threshold}{
	A logical value. If TRUE only markers that fall above the significance threshold alpha.for.pairs (see \code{\link{singlescan}}). If FALSE, all markers are used.
	}
	\item{specific.markers}{
	An optional vector of column numbers or names specifying specific markers to include in the pairscan. 
	}
}
\details{
The markers that were selected for independence can be visualized with the function \code{\link{plotSelectedMarkers}}.
}
\value{
This function returns the data object with two new elements:

\item{geno.for.pairscan}{
	A matrix of the markers to be used in the pair scan. This matrix has the same structure as the original genotype matrix (see \code{\link{read.population}}), but contains only the filtered markers.
	}
\item{covar.for.pairscan}{
	A matrix holding flags for whether each marker should be treated as a covariate. The matrix has one row for each of the filtered markers and one column for each eigentrait being analyzed. The entries of the matrix contain a 1 if the marker is to be used as a covariate for the pair scan and 0 otherwise.
}
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}

\seealso{
\code{\link{singlescan}}, \code{\link{get.linearly.independent}}, \code{\link{plotSelectedMarkers}}
}
\examples{
# not run
# obesity.cross <- select.markers.for.pairscan(obesity.cross, use.pairs.threshold = FALSE)
# plotSelectedMarkers(obesity.cross)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

