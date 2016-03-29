\name{one.singlescan}
\alias{one.singlescan}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Runs one singlescan
}
\description{
This function is called by \code{\link{pairscan.null}} to generate a permuted list of marker effects. The top N markers are selected from this permutation, where N is the number of markers selected for the pairscan. These markers are then permuted in pairs to generate a null distribution for \code{\link{pairscan}}.
}
\usage{
one.singlescan(phenotype.vector, genotype.mat, 
covar.vector = NULL, n.cores = NULL)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{phenotype.vector}{
	A numeric vector holding the phenotype value, or permuted phenotype value, for each individual. 
}
  \item{genotype.mat}{
	A matrix of numeric phenotypes, which are to be associated with the phenotype.
}
  \item{covar.vector}{
	A vector of 0's and 1's, indicating which markers are to be used as covariates in the model.
	}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing. If NULL, the choice is made automatically.
	}
}
\value{
This function returns a table holding the beta, standard error, t statistic and p value for each marker associated with the phenotype. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{pairscan}}, \code{\link{pairscan.null}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

