\name{get.linearly.independent}
\alias{get.linearly.independent}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Find all markers in the genotype matrix that are linearly independent
}
\description{
This function takes in a matrix and tests all columns to verify that they are all linearly independent. If any column is found to be linearly dependent on another, this column is removed from the matrix and its column name is recorded.
}
\usage{
get.linearly.independent(geno.matrix)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{geno.matrix}{
	A matrix contaning genotype information.
	}
}
\details{
This function is called from inside \code{\link{select.markers.for.pairscan}}. It uses correlation between marker pairs to determine which markers are linearly independent. In the restricted space of genotype vectors, two markers will only be linearly non-independent if they have perfect correlation. This function finds the correlations between all pairs of makers. If any correlations are perfect, the first of the markers is picked to discard.

This step is necessary because in the pair scan multiple markers are put into each regression model: The two markers being tested for main effects and interaction, along with any markers being used as covariates. If any markers are linearly dependentent, they do not add any information to the model and R discards them automatically. This filtering step avoids testing extra markers and ensures that the models returned do not have missing values.
}
\value{
This function returns a list with two elements:
\item{independent.matrix}{This is the genotype matrix in which all markers are linearly independent.}
\item{rejected.columns}{A vector of marker names indicating which were removed from the genotype matrix because of linear non-independence.}
}
\references{
Broman, K. W., Wu, H., Sen, S., and Churchill, G. A. (2003). R/qtl: QTL mapping in experimental crosses. Bioinformatics, 19:889-890.
Horn, R. A. and Johnson, C. R. (1985). Matrix Analysis. Cambridge University Press, Cambridge, UK.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{select.markers.for.pairscan}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{algebra}% __ONLY ONE__ keyword per line
\keyword{internal}
