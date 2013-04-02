\name{delete.pheno}
\alias{delete.pheno}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Remove phenotypes from the phenotype matrix
}
\description{
This function deletes an unwanted phenotype or phenotypes from the phenotype matrix.
}
\usage{
delete.pheno(data.obj, phenotypes)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{phenotypes}{
	A vector of either column numbers or column names designating which phenotype or phenotypes should be deleted.
}
}
\value{
This function returns the data object with the specified phenotypes removed.
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
\examples{
data(obesity.cross)
str(obesity.cross)
obesity.cross <- delete.pheno(obesity.cross, "insulin")
str(obesity.cross)
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}
