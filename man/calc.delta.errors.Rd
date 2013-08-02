\name{calc.delta.errors}
\alias{calc.delta.errors}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Error propagation
}
\description{
This function performs error propagation on coefficients and standard errors.
}
\usage{
calc.delta.errors(markers,beta.m,se,beta.cov)
}

\arguments{
  \item{markers}{
	The marker names being tested
	}
\item{beta.m}{
	The main-effects coefficient matrix for the pairwise regression of the given pair.
	}
\item{se}{
	The standard errors for the marker pair.
	}
\item{beta.cov}{
	The model covariance matrix from the pairwise regression
	}
}
\value{
Returns the error propagated coefficients and standard errors for m12 and m21
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}% __ONLY ONE__ keyword per line
