\name{sortCross}
\alias{sortCross}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Sort the genetic markers in the data.obj.
}
\description{
This function sorts genetic markers based on chromosome and chromosomal position.
}
\usage{
sortCross(data.obj, geno.obj)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
	\item{geno.obj}{
	The object in which the genotype matrix and marker information are stored. See \code{\link{read.geno}}.
	}
}
\value{
This function returns the data.obj with the genetic markers sorted by chromosome and position.
}

% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}

