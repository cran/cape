\name{genome.wide.threshold.1D}
\alias{genome.wide.threshold.1D}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Calculate a genome-wide significance threshold for the single-variant scan
}
\description{
This function uses permtuation testing to calculate a genome-wide significance threshold for the single-variant scan. Two user-defined thresholds are calculated: the lower threshold (higher alpha) can be used to determine which variants will be used in the pairwise scan. The higher threshold (lower alpha) is used to determine which variants are used as covariates in the pairwise scan. In each permutation, the phenotype or eigentrait is shuffled and all markers are retested with the permuted phenotype. The regression coefficients are collected from each permutation and the extreme value distribution is used to determine thresholds for the user-defined alpha values.
}
\usage{
genome.wide.threshold.1D(data.obj, n.perm = 1000, alpha.for.pairs = 0.05, alpha.for.covar = 0.01, scan.what = c("eigentraits", "raw.traits"), verbose = FALSE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
  \item{n.perm}{
	The number of permutations to perform. The default is 1000.
}
  \item{alpha.for.pairs}{
	The alpha value used to calculate the threshold for which variants will be used in the pair scan. In \code{pairscan} the use of this threshold is optional.
}
  \item{alpha.for.covar}{
	The alpha value used to calculate the threshold for which variants will be used as covariates in the pairwise scan.
}
  \item{scan.what}{
	A character value that uniquely specifies whether the eigentraits or phenotypes should be scanned.
}
  \item{verbose}{
	A logical value specifying whether the progress of the scan should be printed to the screen. A value of TRUE causes the progress to be printed to the screen.
}
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{link{singlescan}}
\code{\link{pairscan}}
}
\examples{
# not run
# data.obj <- genome.wide.threshold.1D(data.obj)
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

