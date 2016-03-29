\name{impute.missing.geno}
\alias{impute.missing.geno}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Impute missing genotypes in measured markers.
}
\description{
This function imputes missing genotype data using weighted k nearest neighbors imputation.
}
\usage{
impute.missing.geno(data.obj, geno.obj = NULL, 
impute.full.genome = FALSE, k = 10, ind.missing.thresh = 0, 
marker.missing.thresh = 0, 
prioritize = c("ind", "marker", "fewer"), 
max.region.size = NULL, min.region.size = NULL, 
run.parallel = TRUE, verbose = FALSE, n.cores = NULL)
}

\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
\item{geno.obj}{
	The object in which genotype data are stored. See \code{\link{read.geno}}.
	}
\item{impute.full.genome}{
	In CAPE it is possible to scan a subset of the full genotype matrix. This argument indicates whether imputation should be done only on the subset of markers being scanned or on the entire genotype matrix.
	}
\item{k}{
	The number of neighbors to use in k-nearest neighbors imputation.
	}
\item{ind.missing.thresh}{
	A percentage. After imputation the number of data points still missing is assessed. Markers for which large numbers of individuals still have missing data will be removed. This arguement defines the percentage of missing individuals above which a marker will be removed.
	}
\item{marker.missing.thresh}{
	A percentage. After imputation the number of data points still missing is assessed. Individuals for which large numbers of markers still have missing data will be removed. This arguement defines the percentage of missing markers above which an individual will be removed.
	}
\item{prioritize}{
	Markers and individuals with excessive missing data after imputation will be removed. If there is one marker with very low genotyping coverage, it is preferrable to remove that marker rather than all the individuals who are missing a genotype for that marker. This argument allows prioritization of the removal of markers or individuals. The default value is "fewer," meaning that priority goes to whichever dimension (individuals or markers) results in removal of fewer elements.
	}
\item{max.region.size}{
	A number indicating the maximum number of markers to consider in each imputation step. This value defaults to the maximum number of markers on one chromosome.
	}
\item{min.region.size}{
	A number indicating the minimum number of markers to consider in each imputation step. This value defaults to the minimum number of markers on one chromosome.
	}
\item{run.parallel}{
	A logical value indicating whether this process should be run in parallel.
	}
\item{verbose}{
	A logical value indicating whether the progress of the process should be printed to the screen
	}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing. If NULL, the choice is made automatically.
	}
}
\value{
Because both the data.obj and the geno.obj are manipulated by this function, it returns a list of two elements in which the first is the data.obj and the second is the geno.obj. These objects need to be separated to continue with subsequent functions. If no geno.obj was provided to the function, the result it just the data.obj with the genotype matrix embedded.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\seealso{
\code{\link{read.geno}}
\code{\link{read.pheno}}
\code{\link{make.data.obj}}
\code{\link{read.population}}
}