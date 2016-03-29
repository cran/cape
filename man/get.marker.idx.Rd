\name{get.marker.idx}
\alias{get.marker.idx}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Get the column index of markers in the genotype matrix
}
\description{
This function returns the column index in the genotype matrix for each marker specified by the user.
}
\usage{
get.marker.idx(data.obj, markers)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{markers}{
	A vector of either marker names or marker numbers for which the column indices are desired
	}
}
\value{
Returns a vector of numeric values specifying the index of each marker in the genotype matrix.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\seealso{
	\code{\link{get.marker.chr}}
	\code{\link{get.marker.location}}
	\code{\link{get.marker.name}}
	\code{\link{get.marker.num}}
}