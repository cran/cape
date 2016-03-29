\name{linkage.blocks.network}
\alias{linkage.blocks.network}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Calculates linkage blocks using community detection.
	}
\description{
This function calculates linkage blocks based on correlation between markers on a single chromosome. A correlation matrix is converted to a weighted adjacency matrix, and the fastgreedy community detection algorithm from R/igraph is used to define linkage blocks. Adjacent markers in the same community are grouped into a single linkage block.
}
\usage{
linkage.blocks.network(data.obj, geno.obj = NULL, 
collapse.linked.markers = TRUE, threshold.power = 1, 
plot.blocks = FALSE)
}

\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{geno.obj}{
	The object in which the genotype matrix and marker information are stored. See \code{\link{read.geno}}.
	}
\item{collapse.linked.markers}{
	A logical value indicating whether markers should be collapsed into linkage blocks or treated separately.
	}
\item{threshold.power}{
	A value for soft-thresholding the marker-marker correlation matrix before community detection is run.
	}
\item{plot.blocks}{
	A logical value indicating whether the blocks should be plotted. If TRUE, a pdf titled "Recomb.Images.Genotype.Net.Thresh.X.pdf" is generated, where X indicates the soft-thresholding power.
	}
}
\value{
Returns the data object with a new element called linkage.blocks.collapsed if collapse.linked.markers is TRUE, and linkage.blocks.full if collapse.linked.markers is FALSE.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}
\seealso{
\code{\link{get.network}}
}
\references{
Csardi G, Nepusz T: The igraph software package for complex network research, InterJournal, Complex Systems 1695. 2006. http://igraph.org
}