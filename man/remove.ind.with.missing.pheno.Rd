\name{remove.ind.with.missing.pheno}
\alias{remove.ind.with.missing.pheno}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Remove individuals with missing phenotypes
}
\description{
This function removes individuals from a data set if they have missing phenotype data. This is required before the singular value decomposition (SVD) can be performed on the phenotype matrix to obtain eigentraits.
}
\usage{
remove.ind.with.missing.pheno(data.obj)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
}
\value{
The data object is returned after the individuals with missing phenotype data have been removed. A message is also printed to the screen reporting how many individuals were removed.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{get.eigentraits}}
}
\examples{
data(obesity.cross)
obesity.cross <- remove.ind.with.missing.pheno(obesity.cross)
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}
\keyword{manip}% __ONLY ONE__ keyword per line
