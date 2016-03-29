\name{plotSinglescan.heat}
\alias{plotSinglescan.heat}
%- Also NEED an \alias for EACH other topic documented here.
\title{ Plot the results of singlescan as a heatmap
}
\description{
This function plots the results obtained from the single-marker regression performed by singlescan. The effects (\eqn{\beta}) of each regression on each phenotype or eigentrait are plotted as heatmap. This method allows better comparison of effects between phenotypes than \code{\link{plotSinglescan}}.
}
\usage{
plotSinglescan.heat(data.obj, singlescan.obj, 
standardized = TRUE, show.marker.labels = FALSE, 
show.chr.boundaries = TRUE, label.chr = TRUE, 
threshold.above = NULL, color = "blue", 
light.dark = "light", scale.fun = NULL)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
\item{singlescan.obj}{
	The object in which the results of \code{\link{singlescan}} are stored.
	}
	\item{standardized}{
		A logical value. If TRUE, the absolute value of the regression t statistics (\eqn{\beta/\sigma}) are plotted. If FALSE, the raw regression coefficients (\eqn{\beta}) are plotted.
	}
	\item{show.marker.labels}{
		A logical value indicating whether marker names should be printed along the plot axes.
		}
	\item{show.chr.boundaries}{
		A logical value. If TRUE the boundaries of the chromosomes are marked with vertical lines.
	}
	\item{label.chr}{
		A logical value. If TRUE, the chromosomes are labeled with their numbers
	}
	\item{threshold.above}{
		A numeric value. Effects above this value are colored white. This is useful in preventing large effects from swamping out differences between smaller effects.
	}
	\item{color}{
		The color to be used in the raster image.
		}
		\item{light.dark}{
			Either "light", or "dark" to indicate the range of colors that should be plotted.
		}
	\item{scale.fun}{
		The name of function, such as \code{sqrt} or \code{log10}, that will scale the values in the matrix to increase contrast in the image.
	}
}
\value{
Nothing is returned from this function. It produces a heatmap of the effects 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{singlescan}}, \code{\link{plotSinglescan}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{hplot}

