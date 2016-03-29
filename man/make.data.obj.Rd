\name{make.data.obj}
\alias{make.data.obj}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Generate data.obj from pheno.obj and geno.obj
	}
\description{
This function combines information from a phenotype matrix and a genotype object to create a data object in which the results of the analysis will be stored.
}
\usage{
make.data.obj(pheno.obj, geno.obj)
}

\arguments{
  \item{pheno.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{geno.obj}{
	The object in which the genotype matrix and marker information are stored. See \code{\link{read.geno}}.
	}
}
\value{
Returns the data object. Marker information, but not the genotype matrix itself are transferred to the pheno.obj to generate a data.obj. The geno.obj is unchanged and is not returned.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}
\seealso{
\code{\link{read.geno}}
\code{\link{read.pheno}}
\code{\link{read.population}}
}
