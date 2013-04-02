\name{consec.pairs}
\alias{consec.pairs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate a matrix of consecutive pairs
}
\description{
Generates a two-column matrix consisting of consecutive pairs from the input vector. This function in useful for functions that need to be applied to consecutive elements in a vector.
}
\usage{
consec.pairs(elements)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{elements}{
A vector of elements for which consecutive pairs are desired. Elements can be of any class.
}
}
\value{
A two-column matrix is returned. Each row of the matrix consists of two consecutive elements of the input vector. All consecutive pairs are included.
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
\examples{
V <- 1:5
consec.pairs(V)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}
