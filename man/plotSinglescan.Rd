\name{plotSinglescan}
\alias{plotSinglescan}
%- Also NEED an \alias for EACH other topic documented here.
\title{ Plot the results of singlescan
}
\description{
This function plots the results obtained from the single-marker regression performed by singlescan. The effects (\eqn{\beta}) of each regression on each phenotype or eigentrait are plotted as a vertical line. Chromosomes and traits for plotting can be specified.
}
\usage{
plotSinglescan(data.obj, singlescan.obj, 
chr = NULL, traits = NULL, show.alpha.values = NULL, 
standardized = TRUE, show.marker.labels = FALSE, 
mark.covar = FALSE, mark.chr = TRUE, plot.type = "h", 
overlay = FALSE, trait.colors = NULL, 
show.rejected.markers = FALSE, show.selected.markers = FALSE, 
relative.spacing = FALSE, gap.grades = 10, 
min.gap = 1, max.gap = 10, min.gap.genome = 1000, 
max.gap.genome = 1000000, fixed.scale = FALSE, 
ymax = NULL, cex.axis = 0.5, cex.points = 1, 
cex.alpha = 0.5, cex.labels = 2, cex.legend = 0.7)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
	\item{singlescan.obj}{
		The object in which the results from \code{\link{singlescan}} are stored.
	}
	\item{chr}{
		An optional vector indicating which chromosomes should be plotted. If NULL, the default, all chromosomes are plotted.
	}
	\item{traits}{
		An optional vector indicating which traits should be plotted. If NULL, the default, all traits are plotted.
	}
	\item{show.alpha.values}{
		A numeric vector indicating which alpha values specified in \code{\link{singlescan}} should be plotted. If NULL, all alpha values calculated in \code{\link{singlescan}} are plotted.
	}
	\item{standardized}{
		A logical value. If TRUE, the absolute value of the regression t statistics (\eqn{\beta/\sigma}) are plotted. If FALSE, the raw regression coefficients (\eqn{\beta}) are plotted.
	}
	\item{show.marker.labels}{
		A logical value indicating whether marker names should be printed along the plot axes.
		}
	\item{mark.covar}{
		A logical value. If TRUE, the covarates for the pair scan are marked in red. If FALSE, all markers are plotted in black.
	}
	\item{mark.chr}{
		A logical value. If TRUE, alternating chromosomes are shaded in gray to aid in visualizing chromosome boundaries.
		}
	\item{plot.type}{
		A character indicating whether the marker effects should be plotted as vertical lines ("h"), points ("p"), lines ("l"), or both lines and points ("b"). 
	}
	\item{overlay}{
		A logical value indicating whether plots of scans should be overlayed on top of each other or plotted in separate panels.
	}
		\item{trait.colors}{
			If overlay is TRUE, this argument specifies different colors for the ovelayed scans. There are four default colors (black, blue, purple, darkgreen). Alternate or additional colors can be specified manally.
		}
		\item{show.rejected.markers}{
			A logical value. If \code{\link{select.markers.for.pairscan}} has been run, setting this value to TRUE places indicators above markers that have been rejected from being included in the pair scan. show.selected.markers and show.rejected markers cannot be set to TRUE simultaneously.
			}
		\item{show.selected.markers}{
			A logical value. If \code{\link{select.markers.for.pairscan}} has been run, setting this value to TRUE places indicators above markers that have been selected for the pair scan.
			}
	\item{relative.spacing}{
		A logical value indicating whether markers should be plotted at even intervals (FALSE) or spaced relative to their actual genomic locations (TRUE).
	} 
	\item{min.gap}{
		The minimum plotted gap between markers.
		}
	\item{max.gap}{
		The maximum plotted gap between markers.
		}
	\item{min.gap.genome}{
		The minimum gap between markers in genomic distance.
		}
	\item{max.gap.genome}{
		The maximum gap between markers in genomic distance.
		}
	\item{gap.grades}{
		An integer indicating how finely to bin between minimum gap distance and maximum gap distance.
		}
	\item{fixed.scale}{
		A logical value to indicate whether the y axes of multiple plots should have the same minimum and maximum.
		}
	\item{ymax}{
		An optional numeric value indicating a fixed maximum for the y axes.
	}
	\item{cex.axis}{
		A numeric value indicating the size of the font for the axes.
	} 
	\item{cex.points}{
		A numeric value indicating the pointsize or line width for data plotting.
		}
	\item{cex.alpha}{
		A numeric value indicating the font size for writing alpha values.
	}
	\item{cex.labels}{
		A numeric value indicating the font size for labels.
	}
	\item{cex.legend}{
	A numeric value indicating the size of the legend.		
	}

}
\value{
Nothing is returned from this function. It produces a plot of the effects 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{singlescan}}, \code{\link{select.markers.for.pairscan}}
}
\examples{
# plot all markers and both eigentraits
\dontrun{plotSinglescan(obesity.cross)}
# plot only results from chromosomes 1 through 4
\dontrun{plotSinglescan(obesity.cross, chr = c(1:4))}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{hplot}

