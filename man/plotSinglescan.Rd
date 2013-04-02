\name{plotSinglescan}
\alias{plotSinglescan}
%- Also NEED an \alias for EACH other topic documented here.
\title{ Plot the results of singlescan
}
\description{
This function plots the results obtained from the single-marker regression performed by singlescan. The effects (\eqn{\beta}) of each regression on each phenotype or eigentrait are plotted as a vertical line. Chromosomes and traits for plotting can be specified.
}
\usage{
plotSinglescan(data.obj, chr = NULL, traits = NULL, standardized = TRUE, 
mark.covar = TRUE, mark.chr = FALSE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
	\item{chr}{
		An optional vector indicating which chromosomes should be plotted. If NULL, the default, all chromosomes are plotted.
	}
	\item{traits}{
		An optional vector indicating which traits should be plotted. If NULL, the default, all traits are plotted.
	}
	\item{standardized}{
		A logical value. If TRUE, the absolute value of the regression t statistics (\eqn{\beta/\sigma}) are plotted. If FALSE, the raw regression coefficients (\eqn{\beta}) are plotted.
	}
	\item{mark.covar}{
		A logical value. If TRUE, the covarates for the pair scan are marked in red. If FALSE, all markers are plotted in black.
	}
	\item{mark.chr}{
		A lovical value. If TRUE, alternating chromosomes are shaded in gray to aid in visualizing chromosome boundaries.
		}
}
\value{
Nothing is returned from this function. It produces a plot of the effects 
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{singlescan}}
}
\examples{
# not run 
# plot all markers and both eigentraits
# plotSinglescan(obesity.cross)
# plot only results from chromosomes 1 through 4
# plotSinglescan(obesity.cross, chr = c(1:4))
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{hplot}

