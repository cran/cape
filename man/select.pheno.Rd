\name{select.pheno}
\alias{select.pheno}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Select phenotypes for analysis
}
\description{
This function is used to pare down a matrix of phenotypes to only those that will be included in the analysis.
}
\usage{
select.pheno(data.obj, phenotypes)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{phenotypes}{
A vector of either column numbers of the desired phenotypes, or character strings that uniquely identify the columns containing the desired phenotypes.
}
}
\value{
This function returns the data object with a phenotype matrix that contains only the specified phenotypes.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
data(obesity.cross)
obesity.cross <- select.pheno(obesity.cross, c("glucose", "body_weight"))
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}

