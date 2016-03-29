\name{get.marker.num}
\alias{get.marker.num}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Get marker numbers from marker names
}
\description{
This function returns the marker number of each marker specified by name
}
\usage{
get.marker.num(data.obj, markers)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{markers}{
	A vector of marker numbers or names for which names are desired. If the vector contains marker numbers, the marker numbers will be returned unchanged. 
	}
}
\value{
Returns a numeric vector of specifying the number of each specified marker. Markers can either be named with a character string or a number in the genotype matrix. This function allows conversion between them.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\seealso{
	\code{\link{get.marker.chr}}
	\code{\link{get.marker.idx}}
	\code{\link{get.marker.location}}
	\code{\link{get.marker.name}}
}