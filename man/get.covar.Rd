\name{get.covar}
\alias{get.covar}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Get information about covariates
}
\description{
This function returns information about all covariates that originated from either genotype or phenotype.
}
\usage{
get.covar(data.obj, covar = NULL)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{covar}{
		A character string indicating which covariates to return information about.
	}
}
\value{
This function returns a list with the following elements:
\itemize{
\item{covar.names: A vector of covariate names}
\item{covar.type: A vector containing "p" and "g" values to indicate whether each covariate originate as a genetic marker or a phenotype.}
\item{covar.loc: A numeric vector indicating the chromosomal coordinate of each covariate. Phenotypic covariates are given dummy coordinates of 1:n.}
\item{covar.table: A matrix giving the values of each covariate for each individual.}
}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
