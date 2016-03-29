\name{read.pheno}
\alias{read.pheno}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Read in and format data for analysis by \code{\link{cape}}
}
\description{
This function reads in data for \code{cape} analysis and formats it into an object used by other functions in \code{cape}. A single comma-separated file containing both phenotype and genotype data is required. Chromosome and marker locations are required for each marker, and markers are assumed to be in order.
}
\usage{
read.pheno(file.format = c("cape", "csv"), filename = NULL, 
pheno.col = NULL, id.col = 1, delim = ",", na.strings = "-")
}
%- maybe also usage for other objects documented here.
\arguments{
	\item{file.format}{
		A character string indicating which of the accepted formats describes the file to be read in. See Details for specifics.
		}
  \item{filename}{
An optional character string with path name specifying the file to be read in. Omission of this argument will prompt a dialog box for selecting a file.
}
  \item{pheno.col}{
	An optional numeric vector specifying which columns the phenotypes of interest are in. If omitted, all phenotypes are read in.
	}
	\item{id.col}{
		An integer indicating in which column the individual IDs are stored.
		}
\item{delim}{
		A character string indicating the delimeter in the data file. The default indicates a comma-separated file (",").
	}
	\item{na.strings}{
	The symbol used to denote missing data in the file. 	Misspecifying this character can lead to errors in processing the file in which cape misstakenly thinks some phenotypes have character values in them.
	}
}
\details{
Phenotype data can be contained in one of two file formats: cape, or csv. For a description of the cape format, see \code{\link{read.population}}.
The csv format must contain the following:
\itemize{
\item \code{header}: A header labeling each column is required
\item \code{phenotypes}: All phenotypes are required to be numeric. Phenotypes that are not numeric must be coded numerically. For example sex can be coded as [0,1]. Missing values are indicated with the symbol specified by na.strings. The default symbol for na.strings is '-'
}
}
\value{
This function is typically used in conjunction with \code{\link{read.geno}} and \code{\link{make.data.obj}} when genotype data are too large to include in the data.obj, but can be used for small data as well. \code{\link{read.pheno}} initializes the data.obj with the phenotype matrix element "pheno" (see below). After \code{\link{read.geno}} is run, to read in the genotype data separately, \code{\link{make.data.obj}} is run to transfer information from the geno.obj to the initialized data.obj.

\item{pheno}{A matrix containing the phenotype data for the population. Each phenotype is stored in a column, and individuals are stored in rows.}
}

\seealso{
\code{\link{read.population}}, \code{\link{read.geno}}, \code{make.data.obj}
}

% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{IO}

