\name{get.marker.val}
\alias{get.marker.val}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Get marker values
}
\description{
This function returns genotypes for each individual at a given marker. This function is relatively time-consumig and is not recommended for high-throughput use. 
}
\usage{
get.marker.val(data.obj, geno.obj = NULL, markers)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{geno.obj}{
	The object in which the genotype matrix and marker information are stored. See \code{\link{read.geno}}.
	}
\item{markers}{
	A vector of marker numbers or names for which names are desired.	
	}
}
\value{
Returns a matrix listing the individual genotypes of the markers given in the input. Markers can be names either  with a character string or a number.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\seealso{
	\code{\link{get.marker.chr}}
	\code{\link{get.marker.idx}}
	\code{\link{get.marker.location}}
	\code{\link{get.marker.name}}
	\code{\link{get.marker.num}}
}