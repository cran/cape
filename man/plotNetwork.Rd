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
trait = NULL, chr = NULL, show.effect.size = FALSE)
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
  \item{chr}{
	A character vector indicating which chromosomes should be plotted. If NULL, all chromosomes are plotted.
	}
\item{show.effect.size}{
	A logical value indicating whether the height of chromosome main effect bars should vary with the effect size.
	}
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
\seealso{
\code{\link{get.network}}, \code{\link{plotVariantInfluences}}
}
\examples{
# not run
# plotNetwork(obesity.cross, collapsed.net = TRUE)
# plotNetwork(obesity.cross, collapsed.net = FALSE)
# plotNetwork(obesity.cross, collapsed.net = TRUE, trait = "glucose")
# plotNetwork(obesity.cross, collapsed.net = TRUE, chr = c(3:7))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

