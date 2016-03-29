\name{get.marker.location}
\alias{get.marker.location}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Get the chromosomal coordinate of markers
}
\description{
This function returns the chromosomal coordinate of markers specified by the user.
}
\usage{
get.marker.location(data.obj, markers)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{markers}{
	A vector of either marker names or marker numbers for which the chromosomal coordinates are desired.
	}
}
\value{
Returns a vector of numeric values specifying the chromosomal coordinate of each marker in the input vector.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\seealso{
	\code{\link{get.marker.chr}}
	\code{\link{get.marker.idx}}
	\code{\link{get.marker.name}}
	\code{\link{get.marker.num}}
}