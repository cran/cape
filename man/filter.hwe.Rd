\name{filter.hwe}
\alias{filter.hwe}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Filter markers by Hardy-Weinberg equilibrium
}
\description{
This function tests markers for Hardy-Weinberg equilibrium (HWE) and removes markers that are significantly out of HWE.
}
\usage{
filter.hwe(data.obj, geno.obj, p.thresh = 1e-6, run.parallel = TRUE, n.cores = NULL)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{geno.obj}{
	The object in which the genotype object is stored.
	}
  \item{p.thresh}{
	The p value below which markers are considered out of HWE. 
	}
  \item{run.parallel}{
	A logical value to indicate whether to run the process in parallel.
	}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing. If NULL, the choice is made automatically.
	}
}
\value{
This function returns the data object with the marker information pared down to only those markers in HWE.
}
\references{
Add reference for HWE. Hartl and Clark?
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}
\keyword{univar}

