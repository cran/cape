\name{read.population}
\alias{read.population}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Read in and format data for analysis by \code{\link{cape}}
}
\description{
This function reads in data for \code{cape} analysis and formats it into an object used by other functions in \code{cape}. A single comma-separated file containing both phenotype and genotype data is required. Chromosome and marker locations are required for each marker, and markers are assumed to be in order.
}
\usage{
read.population(filename = NULL, pheno.col = NULL, 
geno.col = NULL, delim = ",", na.strings = "-",
check.chr.order = TRUE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{filename}{
An optional character string with path name specifying the file to be read in. Omission of this argument will prompt a dialog box for selecting a file.
}
  \item{pheno.col}{
	An optional numeric vector specifying which columns the phenotypes of interest are in. If omitted, all phenotypes are read in.
}
	\item{geno.col}{
	An optional numeric vector specifying which columns the genotypes of interest are in. If omitted, all genotypes are read in.
	}
	\item{delim}{
		A character string indicating the delimeter in the data file. The default indicates a comma-separated file (",").
	}
	\item{na.strings}{
	The symbol used to denote missing data in the file. 	
	}
	\item{check.chr.order}{
	A logical value indicating whether the order of the chromosomes should be checked. In general, chromosomes should be entered in increasing numerical value. CAPE does not sort chromosomes, and they will be plotted in the order in which they are entered. If the chromosomes have non-numeric and non-X or Y names, and cannot be checked appropriately, or an alternate order is desired, set check.chr.order to FALSE.
	}
}
\details{
All phenotype and genotype data must be contained in a single comma-separated file. The phenotypes should be listed in columns at the beginning of the file, followed by the genotype data. Each row of the file corresponds to one individual. The file must contain the following attributes:
\itemize{
\item \code{header}: A header labeling each column is required
\item \code{chromosomes}: The second line of the file must contain the chromosome on which each marker is found. This line should begin with empty spaces in the phenotype columns followed by a chromosome label for each marker.
\item \code{marker location}: The third line of the file must contain the chromosomal locations of the markers. Like the line of chromosome labels, this line should begin with empty spaces in the phenotype columns followed by a chromosomal position for each marker.
\item \code{phenotypes}: The phenotypes must be listed in the first columns of the file. All phenotypes are required to be numeric. Phenotypes that are not numeric must be coded numerically. For example sex can be coded as [0,1]. Missing values are indicated with the symbol specified by na.strings. The default symbol for na.strings is '-'
\item \code{genotypes}: Genotypes may be coded in one of three different formats: (1) As letters, for example A,H,B, indicating homozygous for allele 1, heterozygous, and homozygous for allele 2 respectively. "H" must be used for heterozygotes, but the other genotypes may be coded with any other letters. (2) As the numbers 0,1,2 indicating homozygous for allele 1, heterozygous, and homozygous for allele 2 respectively. (3) As continuous probabilities of the presence of the reference allele. An individual homozygous for allele 1 would be coded as 0, a heterozygous individual as 0.5, and an individual homozygous for allele 2 as 1. The continuous probabilities allow for uncertainty in genotyping that is not automatically available in the A,H,B or 0,1,2 encodings.
}
}
\value{
The file is converted to a list object that is used as the main argument in most functions. This object is always referred to as data.obj. All data and analysis results will eventually be stored in this object. Upon creation the data.obj contains four elements:
\item{pheno}{A matrix containing the phenotype data for the population. Each phenotype is stored in a column, and individuals are stored in rows.}
\item{geno}{A matrix containing the genotype data for the population. Each genotype is stored in a column, and individuals are stored in rows. Regardless of original format, the genotypes are converted to probabilities for in the data object. Genotypes originally coded as A,H,B for example, will be encoded as 0,0.5,1 respectively.}
\item{chromosome}{A vector containing the chromosome on which each marker is found.}
\item{marker.location}{A vector containing the chromosomal position of each marker.}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{IO}

