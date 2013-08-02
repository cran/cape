\name{singlescan}
\alias{singlescan}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Run the single-variant regression for all phenotypes
}
\description{
This function runs the single-variant regression for either raw phenotypes or eigentraits. It also performs permutation testing and adds two significance thresholds to the data object. One threshold optionally determines the significance cutoff for pairwise testing, and the other optionally determines the significance cutoff for which markers will be used as covariates in the pairwise scan.
}
\usage{
singlescan(data.obj, n.perm = NULL, covar = NULL, 
scan.what = c("eigentraits", "raw.traits"), 
alpha = c(0.01, 0.05), verbose = FALSE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{n.perm}{
	The number of permutations to be performed.
}
  \item{covar}{
	A vector of marker names to be used as covariates. See \code{\link{create.covar}}.
	}
  \item{scan.what}{
	A character string indicating uniquely whether raw traits or eigentraits should be tested.
	}
	\item{alpha}{
	A vector of alpha values for which significant standardized effect values will be calculated.
	}
  \item{verbose}{
	A logical value indicating whether the progress of the scan should be printed to the screen.
}
}
\value{
The data object is returned with two additional elements:
\item{singlescan.results}{A list in which each element corresponds to one trait being scanned. Each element contains a table with the results from the single-marker scan. The table includes, in columns, the effect size marker, the standard error of the effect size, the t statistic from the regression, and the p value. There is one row for each marker.}
\item{covar.flags}{A table indicating which variant, in rows, should be used as a covariate in the scan of each phenotype, in columns. A variant is designated as a covariate with a 1.}
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{select.markers.for.pairscan}}, \code{\link{set.covar}}, \code{\link{get.covar}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{regression}

