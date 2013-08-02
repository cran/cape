\name{rotate.mat}
\alias{rotate.mat}
%- Also NEED an \alias for EACH other topic documented here.
\title{ Rotate a matrix 90 degrees clockwise
}
\description{
This function rotates an image such that when it is plotted by \code{\link{image}} it is plotted in the same orientation that the matrix is printed to the screen.
}
\usage{
rotate.mat(mat)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{mat}{
	A matrix to be rotated
}
}
\value{
A matrix that is rotated 90 degrees clockwise from the original matrix is returned 
}
\examples{
mat <- matrix(1:10, ncol = 10, nrow = 10)
rot.mat <- rotate.mat(mat)
rot.mat
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}
\keyword{manip}% __ONLY ONE__ keyword per line
