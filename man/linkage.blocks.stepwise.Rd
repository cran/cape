\name{linkage.blocks.stepwise}
\alias{linkage.blocks.stepwise}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine linkage blocks based on the correlation between marker genotypes.
}
\description{
This function determines which markers are linked based on the correlation between their genotypes. The first marker on a chromosome is the first marker in the first linkage block. This function steps along the markers in the chromosome until the correlation with the first marker drops below the user-defined threshold. This marker then becomes the first marker in the second block, and the process is repeated until all markers on the chromosome are included in a linkage block. This function is called by \code{\link{get.network}}.
}
\usage{
linkage.blocks.stepwise(data.obj, 
collapse.linked.markers = TRUE, 
r.thresh = 0.7, plot.blocks = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{collapse.linked.markers}{
	A logical value. If TRUE, markers that are linked based on the correlation threshold are collapsed into a single locus. If FALSE, each block contains a single marker.
}
	\item{r.thresh}{
		The minimum Pearson correlation value allowed between the first marker in a linkage block and the last marker in a linkage block.
		}
  \item{plot.blocks}{
	A logical value indicating whether the linkage blocks should be plotted. If TRUE, the function generates plots of the correlation between markers on each chromosome and the boundaries of the linkage blocks determined by the function.
}
}

\value{
This function returns a list. Each element represents a single linkage block and contains the markers included in the linkage block.
}

\seealso{
\code{\link{get.network}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

