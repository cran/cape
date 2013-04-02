\name{get.network}
\alias{get.network}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convert the final results to a form plotted by \code{\link{plotNetwork}} and \code{\link{plotCollapsedVarInf}}
}
\description{
This function converts the significant epistatic interactions to a form that can be plotted as a network. This conversion also optionally condenses the network based on linkage between markers. The degree to which the network is condensed is determined by the argument r2.thresh. This value sets the correlation at which two markers are considered linked. 
}
\usage{
get.network(data.obj, p.or.q = 0.05, 
collapse.linked.markers = TRUE, r2.thresh = 0.8, 
standardize = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{p.or.q}{
			A threshold indicating the maximum adjusted p value considered significant. If an fdr method has been used to correct for multiple testing, this value specifies the maximum q value considered significant.

}
  \item{collapse.linked.markers}{
	A logical value. If TRUE markers are combined into linkage blocks based on correlation. If FALSE, each marker is treated as an independent observation.
	}
  \item{r2.thresh}{
	If collapse.linked.markers is TRUE, this numerical value indicates the Pearson r-squared value at which two markers are considered linked. If two markers are correlated according to the threshold, they are collapsed into a single marker.
}
  \item{standardize}{
	A logical value. If FALSE, the interaction terms are stored as marker to marker influence coefficients. If TRUE, the coefficients are standardized by their standard errors.

}
}
\author{
Vivek M. Philip, Anna L. Tyler, Wei Lu, Justin J. Hendrick, and Gregory W. Carter
}
\seealso{
\code{\link{plotNetwork}}, \code{\link{plotCollapsedVarInf}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

