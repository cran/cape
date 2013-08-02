\name{get.pairs.for.pairscan}
\alias{get.pairs.for.pairscan}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select marker pairs for pairscan based on filters
}
\description{
This function selects which marker pairs can be tested in the pair scan. Even if all markers are linearly independent, some marker pairs may have insufficient recombination between them to populate all genotype combinations. Marker pairs for which genotype combinations have insufficient numbers of individuals are not tested. This function determines which marker pairs have sufficient representation in all genotype combinations. 
}
\usage{
get.pairs.for.pairscan(geno, min.per.genotype = NULL,
max.pair.cor = NULL, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{
	A two dimensional genotype matrix with rows containing individuals and columns containing markers. Each entry is a value between 0 and 1 indicating the genotype of each individual at each marker. 
}
  \item{min.per.genotype}{
	The minimum number of individuals allowable per genotype. If for a given marker pair, one of the genotypes is underrepresented, the marker pair is not tested. If this value is NULL, max.pair.cor must have a numeric value.
}
	\item{max.pair.cor}{
		A numeric value between 0 and 1 indicating the maximum Pearson correlation that two markers are allowed. If the correlation between a pair of markers exceeds this threshold, the pair is not tested. If this value is set to NULL, min.per.genotype must have a numeric value.
	}
  \item{verbose}{
	A logical value. If TRUE, the script prints a message to the screen to indicate that it is running. If FALSE, no message is printed.
	}
}
\value{
This function returns a two-column matrix of marker pairs. This matrix is then used as an argument in \code{\link{one.pairscan}} to specify which marker pairs should be tested. 
}

\seealso{
\code{\link{one.pairscan}}
}
\examples{
#This example finds marker pairs to test in a randomly
#generated genotype matrix
possible.genotypes <- c(0, 0.5, 1)
genotype.matrix <- matrix(sample(possible.genotypes, 1000, replace = TRUE), nrow = 100, ncol = 10)
colnames(genotype.matrix) <- paste("marker", letters[1:10], sep = "_")
pairs.which <- get.pairs.for.pairscan(genotype.matrix, min.per.genotype = 6, verbose = TRUE)
print(pairs.which)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arith}
\keyword{internal}
