\name{plotPhenoCor}
\alias{plotPhenoCor}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot correlations between phenotype pairs.
}
\description{
This function plots pairs of phenotypes against each other to examine phenotype correlations. CAPE performs best when phenotypes are moderately correlated (0.4 < r < 0.8). This script can help prioritize phenotypes for use in CAPE. Phenotype pairs are plotted in a panel in which the lower triangle shows points colored by the indicated covariate, and the upper triangle shows correlation statistics.
}
\usage{
plotPhenoCor(data.obj, pheno.which = NULL, 
color.by = NULL, group.labels = NULL, 
text.cex = 1, pheno.labels = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{pheno.which}{
	A vector of character strings indicating which phenotypes should be plotted. If pheno.which is NULL, all phenotypes will be plotted.
		}
	\item{color.by}{
	A vector of character strings indicating how the points should be colored. Points are colored based on covariates. These covariates must be specified using \code{\link{pheno2covar}} or \code{\link{marker2covar}} before they can be recognized by this function.
	}
	\item{group.labels}{
		A vector of character strings naming the groups that are specified by the covariates. For example, if coloring by the covariate "sex," coded as 0 for females and 1 for males, the group labels would be c("female", "male"). The vector must be the same length as the number of groups, and must match the numerical order of the groups. 
		}
	\item{text.cex}{
		A numerical value indicating the size of the text showing the correlation values between phenotypes. 
		}
	\item{pheno.labels}{
		An optional vector of character strings giving alternate names for the phenotypes being plotted.
	}
}

\seealso{
\code{\link{pheno2covar}}, \code{\link{marker2covar}}, \code{\link{plotPheno}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

