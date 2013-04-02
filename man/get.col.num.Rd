\name{get.col.num}
\alias{get.col.num}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Use column titles to retrieve column numbers
}
\description{
This function is used internally to tranlate column names to column numbers.
}
\usage{
get.col.num(data.mat, col.which = NULL, warn = TRUE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.mat}{
	A matrix with column names
}
  \item{col.which}{
	The names of the columns for which numbers are wanted
}
  \item{warn}{
	A logical indicating whether the use should be warned if any of the entries in col.which are not found.
}
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
\examples{

mat <- matrix(rnorm(50), ncol = 5)
colnames(mat) <- letters[1:5]
print(mat)
col.num <- get.col.num(mat, c("c", "e"))
print(col.num)

col.num <- get.col.num(mat, c("b", "g"))
print(col.num)

}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

