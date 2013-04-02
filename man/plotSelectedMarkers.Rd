\name{plotSelectedMarkers}
\alias{plotSelectedMarkers}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot markers selected for the pair scan.
}
\description{
This function is for visualizing the results of \code{\link{select.markers.for.pairscan}}. It plots the results from the
single-locus scan in the same type of plot as \code{\link{plotSinglescan}}, and indicates which markers were selected for 
the pair scan based on their linear independence from all other markers in the cross.
}
\usage{
plotSelectedMarkers(data.obj, mark.which = c("selected", "rejected"), 
chr = NULL, traits = NULL, standardized = TRUE, mark.covar = FALSE,
mark.chr = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
  \item{mark.which}{
	A character value indicating whether the markers selected for the pair scan or the markers rejected from the pair scan should be indicated. If set to NULL, no markers are indicated, and the resulting image is the same as in \code{\link{plotSinglescan}}.
	}
  \item{chr}{
	An optional numeric vector indicating which chromosomes should be plotted. The default behavior is to plot all chromosomes.
	}
	\item{traits}{
	An optional vector of trait names indicating which should be plotted. The default behavior is to plot all traits.
	}
	\item{standardized}{
		A logical value. If TRUE, the values plotted are the beta coefficients divided by the standard error (\eqn{\beta/\sigma}). If FALSE, the raw beta coefficients from the regression are plotted.
	}
	\item{mark.covar}{
		A logical value. If TRUE, the markers that are specified as covariates are colored red.
	}
	\item{mark.chr}{
		A lovical value. If TRUE, alternating chromosomes are shaded in gray to aid in visualizing chromosome boundaries.
	}
}
\value{
No values are returned.
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}

\seealso{
\code{\link{select.markers.for.pairscan}}, \code{\link{plotSinglescan}}
}
\examples{
# not run
# obesity.cross <- select.markers.for.pairscan(obesity.cross)
# plotSelectedMarkers(obesity.cross, mark.which = "selected")
# plotSelectedMarkers(obesity.cross, mark.which = "rejected")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

