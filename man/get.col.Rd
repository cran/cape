\name{get.col}
\alias{get.col}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Retrieve colors for use in plotting.
}
\description{
This function is used internally to retrieve colors for plotting.
}
\usage{
get.col(col.name, light.dark = c("light", "dark"))
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{col.name}{
	A character string indicating the color values to be retrieved. One of "green", "purple", "orange", "blue", "brown", "yellow", or "gray".
}
  \item{light.dark}{
	A character string indicating whether the values returned should be in a light or dark range.
	}
}
\value{
Returns a vector of three hex values specifying the colors requested.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

