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
get.network(data.obj, p.or.q = 0.05, min.std.effect = 0,
collapse.linked.markers = TRUE, r.thresh = 0.5, 
verbose = FALSE, plot.linkage.blocks = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{p.or.q}{
			A numerical threshold indicating the maximum adjusted p value considered significant. If an fdr method has been used to correct for multiple testing, this value specifies the maximum q value considered significant.

}
	\item{min.std.effect}{
		A numerical threshold indicating the absolute value of the minimum standard effect size to be shown in the plot. The default value of 0 performs no thresholding.
	}
  \item{collapse.linked.markers}{
	A logical value. If TRUE markers are combined into linkage blocks based on correlation. If FALSE, each marker is treated as an independent observation.
	}
  \item{r.thresh}{
	If collapse.linked.markers is TRUE and linkage.method is set to "genotype", this numerical value indicates the Pearson r value at which two markers are considered linked. If two markers are correlated according to the threshold, they are collapsed into a single marker.
}
	\item{verbose}{
		A logical value indicating whether the function progress should be printed to the screen.
	}
	\item{plot.linkage.blocks}{
		A logical value indicating whether the chromosomes should be plotted with their linkage blocks delineated. The type of plot produced differs depending on which choice is specified by linkage.method.
	}
}
	\details{
		This function calls defineds linkage blocks based on the correlation between adjacent markers. The first marker on the chromosome is the first marker of the initial linkage block. The algorithm steps through markers individually comparing each to the first marker in the block. When a marker has a correlation lower than the specified threshold when compared to the initial marker in the block, a boundary is drawn, and the marker becomes the first marker of the next block. This process is repeated until all markers on the chromosome are assigned to a block.
	}
\seealso{
\code{\link{plotNetwork}}, \code{\link{plotCollapsedVarInf}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

