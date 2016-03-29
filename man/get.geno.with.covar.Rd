\name{get.geno.with.covar}
\alias{get.geno.with.covar}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Retrieve the genotype matrix with the covariates appended.
}
\description{
This function is similar to get.geno, but it returns the genotype matrix with covariates appended for greater ease of testing all marker pairs including covariates. Covariates derived from genetic markers are placed in the appropriate location based on chromosome coordinates.
}
\usage{
get.geno.with.covar(data.obj, geno.obj = NULL, 
g.covar = TRUE, p.covar = TRUE, 
for.pairscan = TRUE)
}

\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{geno.obj}{
	The object in which genotype data are stored. See \code{\link{read.geno}}.
	}
\item{g.covar}{
	A logical value indicating whether covariates derived from genetic markers should be included.
	}
\item{p.covar}{
	A logical value indicating whether covariates derived from phenotypes should be included.
	}
\item{for.pairscan}{
	A logical value indicating if the full genotype matrix should be returned, or only the markers being evaluated in the pairscan.
	}
}
\value{
Returns a matrix of genotype values for the individuals and markers specified in the data object with any covariates included.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\seealso{
\code{\link{get.geno}}
}
\keyword{internal}