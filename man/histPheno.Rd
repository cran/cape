\name{histPheno}
\alias{histPheno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot histograms of phenotypes.
}
\description{
This function plots histograms of phenotypes.
}
\usage{
histPheno(data.obj, pheno.which = NULL, 
pheno.labels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
	\item{pheno.which}{
	A vector of character strings indicating which phenotypes should be plotted. If pheno.which is NULL, all phenotypes will be plotted.
		}
	\item{pheno.labels}{
		An optional vector of character strings giving alternate names for the phenotypes being plotted.
	}
}

\seealso{
\code{\link{plotPhenoCor}}, \code{\link{plotPheno}}, \code{\link{qqPheno}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

