\name{writePopulation}
\alias{writePopulation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write out a cape data object to .csv format.
}
\description{
This function takes in a data object and writes it to a .csv file. This file can be read in by \code{\link{read.population}}.
}
\usage{
writePopulation(data.obj, filename, 
geno.to.convert = NULL, conversion = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{filename}{
	A string specifying the name of the file to be written.
}
	\item{geno.to.convert}{
		A vector of genotypes in the data object that are to be converted to another form. For example, to convert genotypes coded as (0, 0.5, 1), to ("A", "H", "B"), set geno.to.convert to c(0, 0.5, 1) and conversion to c("A", "H", "B").
	}
	\item{conversion}{
		A vector the same length as geno.to.convert. The values in this vector will be written in place of the values in geno.to.convert in the final file.
		}
}
\value{
This function writes a data object to a .csv file that can be read by \code{\link{read.population}}. Nothing is returned.
}
\seealso{
	\code{\link{read.population}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}

