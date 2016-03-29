\name{remove.ind}
\alias{remove.ind}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Remove individuals from the data.obj
}
\description{
This function removes individuals from a data set.
}
\usage{
remove.ind(data.obj, ind.which)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{ind.which}{
		A numerical vector indicating which individuals should be removed.
	}
}
\value{
The data object is returned with the specified individuals removed from the phenotype matrix.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{remove.markers}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}% __ONLY ONE__ keyword per line
