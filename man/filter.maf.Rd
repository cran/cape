\name{filter.maf}
\alias{filter.maf}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Filter markers by minor allele frequency.
}
\description{
This function tests markers for minor allele frequency and removes markers below the user-set threshold.
}
\usage{
filter.maf(data.obj, geno.obj, maf = 0.05)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{geno.obj}{
	The object in which the genotype object is stored.
	}
  \item{maf}{
	The minor allele frequency below which markers are removed.
	}
}
\value{
This function returns the data object with the marker information pared down to only those markers with minor allele frequencies greater than that specified by the user.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}
\keyword{univar}

