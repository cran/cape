\name{select.by.chr}
\alias{select.by.chr}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Subset a cross object to include only specified chromosomes.
}
\description{
This function subsets a data object to include only specified chromosomes.
}
\usage{
select.by.chr(data.obj, chr, include.covariates = TRUE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{chr}{
	A vector of desired chromosomes
}
	\item{include.covariates}{
	Covariates are stored in the genotype matrix. To include these in the final subset, set include.covariates to TRUE. If include.covariates is FALSE, covariates will not automatically be included in the final genotype matrix.
	}
}
\value{
This function returns the data object with the genotype matrix pared down to only the desired chromosomes.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{select.by.ind}}, \code{\link{create.covar}}
}
\examples{
data(obesity.cross)
obesity.cross <- create.covar(obesity.cross, "mom")
# Subset the cross to only include chromosome 6 and the covariates
obesity.cross <- select.by.chr(obesity.cross, c(0, 6))
str(obesity.cross)
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}
