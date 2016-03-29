\name{plotNetwork}
\alias{plotNetwork}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the final epistatic network
}
\description{
This function plots the final results with a using a different layout than \code{\link{plotVariantInfluences}}. Instead of the adjacency matrix, this function plots interactions between markers with arrows indicating the direction of influence. The function \code{\link{get.network}} must be run before plotting the network.
}
\usage{
plotNetwork(data.obj, collapsed.net = TRUE, 
trait = NULL, phenotype.labels = NULL, 
main.lwd = 4, inter.lwd = 3, label.cex = 1.5, 
percent.bend = 15, chr.gap = 1, label.gap = 5, 
positive.col = "brown", negative.col = "blue")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
  \item{collapsed.net}{
	A logical value. If TRUE, a network condensed by linkage (see \code{\link{get.network}}) is plotted. If FALSE, the full network is plotted. 
}
  \item{trait}{
	A character vector indicating which traits should be plotted. If NULL, all traits are plotted.
}
  \item{phenotype.labels}{
	A vector of strings listing alternate names for the phenotypes.
	}
\item{main.lwd}{
	A numeric value specifying the line width of the main effects bars surrounding the chromosome bars.
	}
\item{inter.lwd}{
	A numeric value specifying the line width of the interaction arrows drawn inside the chromosome circle
	}
\item{label.cex}{
	A numberic value specifying the size of the chromomsome labels.
	}
\item{percent.bend}{
	A numeric value between 0 and 100 specifying how much curvature should be added to the arrows within the circle.
	}
\item{chr.gap}{
	A numeric value between 0 and 100 indicating what percentage of the circle should be devoted to individual gaps between chromosomes.
	}
\item{label.gap}{
	A numeric value between 0 and 100 indicating what percentage of the circle should be devoted to the phenotype labels.
	}
\item{positive.col}{
	 The color in which to plot the postive values. Colors can be one of the following: "green", "purple", "orange", "blue", "brown", "yellow", or "gray."

	}
\item{negative.col}{
	 The color in which to plot the negative values. Colors can be one of the following: "green", "purple", "orange", "blue", "brown", "yellow", or "gray."
	}


}
\seealso{
\code{\link{get.network}}, \code{\link{plotVariantInfluences}}
}
\examples{
\dontrun{
	plotNetwork(obesity.cross, collapsed.net = TRUE)
	plotNetwork(obesity.cross, collapsed.net = FALSE)
	plotNetwork(obesity.cross, collapsed.net = TRUE, trait = "glucose")
	}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

