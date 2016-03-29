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
standardized = TRUE, min.std.effect = 0, 
collapse.linked.markers = TRUE, verbose = FALSE, 
plot.linkage.blocks = FALSE, threshold.power = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{p.or.q}{
			A numerical threshold indicating the maximum adjusted p value considered significant. If an fdr method has been used to correct for multiple testing, this value specifies the maximum q value considered significant.

}
	\item{standardized}{
		A logical value indicating whether values placed in weighted adjacency matrix are standardized values or unstandardized values.
	}
	\item{min.std.effect}{
		A numerical threshold indicating the absolute value of the minimum standard effect size to be shown in the plot. The default value of 0 performs no thresholding.
	}
  \item{collapse.linked.markers}{
	A logical value. If TRUE markers are combined into linkage blocks based on correlation. If FALSE, each marker is treated as an independent observation.
	}
	\item{verbose}{
		A logical value indicating whether the function progress should be printed to the screen.
	}
	\item{plot.linkage.blocks}{
		A logical value indicating whether the chromosomes should be plotted with their linkage blocks delineated. The type of plot produced differs depending on which choice is specified by linkage.method.
	}
	\item{threshold.power}{
		The linkage blocks are calculated by finding communities in a correlation matrix of genetic markers. This power indicates a soft-thresholding power for the correlation matrix. The higher the power, the more blocks will be generated.
	}
}
	\details{
		This function delineates linkage blocks based on the correlation between adjacent markers. A correlation matrix is first calculated between all marker pairs on a single chromosome. This similarity matrix is used to construct a weighted network, and the fastgreedy community detection algorithm from the R package igraph is used to detect individual communities within the network. Each community of contiguous markers is designated as a distinct linkage block.
	}
\seealso{
\code{\link{plotNetwork}}, \code{\link{plotCollapsedVarInf}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\references{
Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. http://igraph.org
}
\keyword{manip}

