\name{select.by.chr}
\alias{select.by.chr}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Subset a cross object to include only specified chromosomes.
}
\description{
This function subsets a data object to include only specified chromosomes.
}
\usage{
select.by.chr(data.obj, chr)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
  \item{chr}{
	A vector of desired chromosomes
	}
}
\value{
This function returns the data object with the genotype matrix pared down to only the desired chromosomes.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{select.by.ind}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}
