\name{get.marker.chr}
\alias{get.marker.chr}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Get chromosome assignments for a vector of markers.
}
\description{
This function returns the chromosome that each marker in the user-supplied vector is on.
}
\usage{
get.marker.chr(data.obj, markers)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{markers}{
	A vector of either marker names or marker numbers for which chromosome assignments are desired.
	}
}
\value{
Returns a vector of numeric values specifying which chromosome each marker in the input vector resides on.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\seealso{
	\code{\link{get.marker.idx}}
	\code{\link{get.marker.location}}
	\code{\link{get.marker.name}}
	\code{\link{get.marker.num}}
}