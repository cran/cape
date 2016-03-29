\name{plotCollapsedVarInf}
\alias{plotCollapsedVarInf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot variant-to-variant influences
}
\description{
This function plots the final network as does \code{\link{plotVariantInfluences}}, but it plots the network condensed by linkage blocks as performed by \code{\link{get.network}}. This function can only be run after running \code{\link{get.network}}
}
\usage{
plotCollapsedVarInf(data.obj, expand.labels = FALSE, 
all.markers = FALSE, scale.effects = c("none", "log10", "sqrt"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
	\item{expand.labels}{
		A logical value. If FALSE, the markers are labeled as linkage blocks ("block1", "block2" and so forth). If TRUE, the block labels are expanded to show allmarker names included in each block.
	}
	\item{all.markers}{
		A logical value. If TRUE all markers are plotted. If FALSE only markers tested in the pair scan are plotted.
	}
	\item{scale.effects}{
		A string indicating a scaling function by which the effects should be scaled. This is useful in increasing contrast between effects with large variance.
	}
}

\seealso{
\code{\link{pairscan}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

