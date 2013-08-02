\name{get.covar}
\alias{get.covar}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Use a threshold to automatically select covariates for the pairscan
}
\description{
This function uses a covariate threshold to calculate which markers will be used as covariates in the pair scan.
}
\usage{
get.covar(data.obj, covar.thresh = NULL)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{covar.thresh}{
		A numerical value indicating the standardized effect size (\eqn{\beta/\sigma}) above which a marker is considered to have a strong enough effect to be used as a covariate.
	}
}
\value{
This function operates on the element in data.obj called \code{covar.flags}. This element is a table with one row for each marker and one column for each trait being analyzed. The entry for each marker contains a 1 if it is to be used as a covariate and a 0 otherwise.
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
\seealso{
\code{\link{set.covar}}, \code{\link{singlescan}}, \code{\link{select.markers.for.pairscan}}, \code{\link{pairscan}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{arith}

