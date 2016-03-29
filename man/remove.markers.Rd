\name{remove.markers}
\alias{remove.markers}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Remove markers from the data.obj
}
\description{
This function removes genetic markers from a data set.
}
\usage{
remove.markers(data.obj, markers.which)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{markers.which}{
		A vector indicating the names or ID numbers of which markers should be removed.
		}
}
\value{
The data object is returned with the specified markers removed.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{remove.ind}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}% __ONLY ONE__ keyword per line
