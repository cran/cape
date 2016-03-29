\name{pheno2covar}
\alias{pheno2covar}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Create a covariate from a phenotype.
	}
\description{
Covariates can be derived from either genetic markers or phenotypic values. This function generates covariates from phenotypes. These covariates are treated distinctly from those generated from genetic markers. For example, covariates derived from genetic markers are used in kinship calculations, whereas covariates derived from phenotypes are not. Covariates derived from phenotypes are also plotted after all genetic markers.
}
\usage{
pheno2covar(data.obj, pheno.which)
}

\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{pheno.which}{
	A vector of character strings indicating which phenotypes should be converted to covariates.
	}
}
\value{
Returns the data object with two additional elements.
\item{g.covar.table}{A table of the covariates and their values.}
\item{p.covar}{A vector of the names of the phenotypic covariates}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}
\seealso{
\code{\link{marker2covar}}
}
