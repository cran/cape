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
select.markers.for.pairscan(data.obj, geno.obj = NULL, 
singlescan.obj, alpha.thresh = NULL, 
t.thresh = NULL, specific.markers = NULL, 
num.markers = NULL, step.size = NULL, 
tolerance = 10, verbose = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
	\item{geno.obj}{
		The object in which the genotype data is stored. See \code{\link{read.geno}}.
		}
	\item{singlescan.obj}{
		The object in which results from \code{\link{singlescan}} are stored.
	}
  \item{alpha.thresh}{
	A numeric threshold for selecting markers based on significance at the given alpha value.
	}
	\item{t.thresh}{
	A numerical value indicating a standardized effect size below which markers are rejected for the pair scan.
	}
	\item{specific.markers}{
	An optional vector of column numbers or names specifying specific markers to include in the pairscan. 
	}
	\item{num.markers}{
		An optional number. If specified, the algorithm attempts to find a standardized effect threshold that will give the number of markers specified. This can take a long time, since linear non-independence must be calculated for each thresholded set of markers.
	}
	\item{step.size}{
	If num.markers is specified, this argument determines how the t.thresh is changed when searching for markers to include in the pairscan.
	}
	\item{tolerance}{
		If num.markers is specified, tolerance determines how close the algorithm will get to the number of desired markers. If 100 markers are desired, and tolerance is set to 10, the algorithm will stop if it finds a threshold giving between 90 and 110 markers.
	}
	\item{verbose}{
		A logical value indicating whether the progress of the threshold finding algorithm should be reported.
	}
}
\details{
The markers that were selected for independence can be visualized with the function \code{\link{plotSinglescan}} with either show.selected.markers or show.rejected.markers set to TRUE.
}
\value{
This function returns the data object with an added genotype matrix to be used in the pair scan. This matrix has the same structure as the original genotype matrix, but contains only the filtered markers.
}

\seealso{
\code{\link{singlescan}}, \code{\link{plotSinglescan}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

