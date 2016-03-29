\name{get.pheno}
\alias{get.pheno}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Retrieve the genotype matrix.
}
\description{
This function retrieves the phenotype matrix either from the data object. 
}
\usage{
get.pheno(data.obj, scan.what = c("eigentraits", "raw.traits"))
}

\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{scan.what}{
	A character string defining whether the eigentraits or phenotypes are desired.
	}
}
\value{
Returns a matrix of phenotype values.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\seealso{
\code{\link{read.geno}}
\code{\link{read.pheno}}
\code{\link{make.data.obj}}
\code{\link{read.population}}
\code{\link{get.geno}}
}