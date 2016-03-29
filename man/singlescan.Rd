\name{singlescan}
\alias{singlescan}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Run the single-variant regression for all phenotypes
}
\description{
This function runs the single-variant regression for either raw phenotypes or eigentraits using a kinship correction if desired. It also performs permutation testing and calculates significance thresholds.
}
\usage{
singlescan(data.obj, geno.obj = NULL, n.perm = NULL, 
covar = NULL, alpha = c(0.01, 0.05), 
scan.what = c("eigentraits", "raw.traits"), 
use.kinship = FALSE, kin.full.geno = FALSE, 
run.parallel = TRUE, sample.kinship = FALSE, 
num.kin.samples = 1000, n.per.sample = 10, 
verbose = FALSE, overwrite.alert = TRUE, n.cores = NULL)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
	\item{geno.obj}{
	The object in which the genotype matrix and marker information are stored. See \code{\link{read.geno}}.
	}
  \item{n.perm}{
	The number of permutations to be performed.
	}
  \item{covar}{
	A vector of marker names to be used as covariates. See \code{\link{pheno2covar}} and \code{\link{marker2covar}}.
	}
	\item{alpha}{
	A vector of alpha values for which significant standardized effect values will be calculated.
	}
  \item{scan.what}{
	A character string indicating uniquely whether raw traits or eigentraits should be tested.
	}
	\item{use.kinship}{
		A logical value indicating whether a kinship correction should be implemented.
	}
	\item{kin.full.geno}{
		A logical value indicating whether the entire genotype matrix should be used to calculate kinship corrections. If FALSE, only the subset of genetic markers being scanned will be used to calculate kinship corrections.
	}
	\item{run.parallel}{
		A logial value indicating whether to run calculations in parallel.
		}
	\item{sample.kinship}{
		A logical value indicating whether kinship matrices should be sampled (TRUE) or calculated directly (FALSE). Sampling kinshp matrices can be faster in the case of very large genotype matrices.
	}
	\item{num.kin.samples}{
		If sample.kinship is TRUE, this integer indicates how many samples should be used to calculate kinship matrices.
	}
	\item{n.per.sample}{
		If sample.kinship is TRUE, this integer indicates how many markers should be used in each sample of the genotype matrix for calculating kinship matrices.
	}
  \item{verbose}{
	A logical value indicating whether the progress of the scan should be printed to the screen.
	}
\item{overwrite.alert}{
	If TRUE, this triggers an alert warning the user that the output of this function should be saved separately from the data.obj.
	}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing. If NULL, the choice is made automatically.
	}
}
\value{
IMPORTANT NOTE: Prior versions of this function modified the data.obj. The new version creates a separate singlescan.obj, and the output of this function should NOT be assigned to the data.obj.
The singlescan.obj is a list containing the following four elements:
\item{alpha}{A vector containing the alpha values requested by the user.}
\item{alpha.thresh}{A list containing the same number of elements as alpha values requested. For each value, a significance threshold in terms of t statistic is given.}
\item{singlescan.results}{A list in which each element corresponds to one trait being scanned. Each element contains a table with the results from the single-marker scan. The table includes, in columns, the effect size marker, the standard error of the effect size, the t statistic from the regression, and the p value. There is one row for each marker.}
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{plotSinglescan}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{regression}

