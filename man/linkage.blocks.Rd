\name{linkage.blocks}
\alias{linkage.blocks}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Determine linkage blocks
}
\description{
This function determines which markers are linked based on a user-defined correlation threshold and collapses the markers into linkage blocks. This function is called by \code{\link{get.network}}.
}
\usage{
linkage.blocks(data.obj, p.or.q = 0.05, collapse.linked.markers = TRUE, r2.thresh = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{p.or.q}{
		A threshold indicating the maximum adjusted p value considered significant. If an fdr method has been used to correct for multiple testing, this value specifies the maximum q value considered significant. Only marker pairs with p or q values below this threshold will be plotted.
}
  \item{collapse.linked.markers}{
	A logical value. If TRUE, markers that are linked based on the correlation threshold are collapsed into a single locus. If FALSE, each block contains a single marker.
}
  \item{r2.thresh}{
	The minimum Pearson r-squared value at which two markers are considered to be in the same linkage block.
}
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
\seealso{
\code{\link{get.network}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

