\name{error.prop}
\alias{error.prop}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Estimate Errors of Regression Coefficients
}
\description{
This function uses error propagation formulas for quantities computed from regression coefficients to estimate the error for all regression coefficients.
}
\usage{
error.prop(data.obj, pairscan.obj, perm = FALSE, verbose = FALSE, n.cores = NULL)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{pairscan.obj}{
	The object in which the results from the pairscan are stored. See \code{\link{pairscan}}.
	}
  \item{perm}{
A logical value to indicate whether error propagation should be performed on the test statistics (FALSE) or the permuted test statistics (TRUE).
}
  \item{verbose}{
	A logical value to indicate whether the progress of the function should be printed to the screen.
	}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing. If NULL, the choice is made automatically.
	}
}
\value{
This function returns the data object with a new list element: var.to.var.influences if perm is set to FALSE and var.to.var.influences.perm if perm is set to TRUE. These tables include the errors calculated for the marker1 to marker2 influences as well as the marker2 to marker1 influences. These results are used by \code{\link{calc.p}} to calculate empirical p values.
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010

Bevington PR (1969) Data Reduction and Error Analysis for the Physical Sciences. New York: McGraw-Hill.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{math}
\keyword{symbolmath}

