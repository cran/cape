\name{get.geno}
\alias{get.geno}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Retrieve the genotype matrix.
}
\description{
This function retrieves the genotype matrix either from the data object or the genotype object. If the genotype matrix is held in a separate genotype object, the marker and individual information in the data object is used to filter the genotype matrix to the appropriate markers and individuals.
}
\usage{
get.geno(data.obj, geno.obj)
}

\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{geno.obj}{
	The object in which genotype data are stored. See \code{\link{read.geno}}.
	}
}
\value{
Returns a matrix of genotype values for the individuals and markers specified in the data object.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\seealso{
\code{\link{read.geno}}
\code{\link{read.pheno}}
\code{\link{make.data.obj}}
\code{\link{read.population}}
\code{\link{get.pheno}}
}