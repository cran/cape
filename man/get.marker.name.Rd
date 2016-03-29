\name{get.marker.name}
\alias{get.marker.name}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Get marker names from marker numbers
}
\description{
This function returns the name of markers specified by number
}
\usage{
get.marker.name(data.obj, markers)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{markers}{
	A vector of marker numbers or names for which names are desired. If the vector contains marker names, the marker names will be returned unchanged. 
	}
}
\value{
Returns a vector of character strings specifying the name of each marker designated by a number or name in the input vector. Markers can either be named with a character string or a number in the genotype matrix. This function allows conversion between them.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\seealso{
	\code{\link{get.marker.chr}}
	\code{\link{get.marker.idx}}
	\code{\link{get.marker.location}}
	\code{\link{get.marker.num}}
}