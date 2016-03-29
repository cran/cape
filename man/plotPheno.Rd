\name{plotPheno}
\alias{plotPheno}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot phenotype values by individual.
}
\description{
This function plots phenotype values by individual, and can color points based on covariates if desired.
}
\usage{
plotPheno(data.obj, pheno.which = NULL, 
color.by = NULL, group.labels = NULL, 
pheno.labels = NULL)
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
	\item{pheno.labels}{
		An optional vector of character strings giving alternate names for the phenotypes being plotted.
	}
}

\seealso{
\code{\link{pheno2covar}}, \code{\link{marker2covar}}, \code{\link{plotPhenoCor}}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

