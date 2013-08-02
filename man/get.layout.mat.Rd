\name{get.layout.mat}
\alias{get.layout.mat}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Get the best layout matrix for a given number of panes per page.
}
\description{
This function is for automatically deciding how to lay out multiple images per page. It takes as an argument the number of images the user wants on a single page and returns the most square matrix possible.
}
\usage{
get.layout.mat(num.panes, type = c("landscape", "upright"))
}

\arguments{
  \item{num.panes}{
	The number of images to plot per page.
	}
\item{type}{
	A character vector specifying whether the layout matrix should have more columns than rows ("landscape"), or more rows than columns ("upright"), if the resulting matrix is not square.
	}
}
\value{
Returns a matrix that can be used as a layout matrix by the function \code{\link{layout}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}% __ONLY ONE__ keyword per line
