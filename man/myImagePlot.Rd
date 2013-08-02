\name{myImagePlot}
\alias{myImagePlot}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Used in plotting results of pair scan.
}
\description{
This is an internal plotting function that plots matrix images on a yellow and blue scale with a color bar. It is called by \code{\link{plotPairscan}}
}
\usage{
myImagePlot(x, ...)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{x}{
	The matrix to be plotted as an image
}
  \item{\dots}{
	Additional plotting arguments, such as x and y labels. MyImagePlot also looks for min.x and max.x arguments to specify the scale of z values on which the x matrix should be plotted.
}
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{plotPairscan}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}
\keyword{hplot}% __ONLY ONE__ keyword per line
