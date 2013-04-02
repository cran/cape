\name{create.covar}
\alias{create.covar}
%- Also NEED an \alias for EACH other topic documented here.
\title{
	Designate a phenotype or treatment as a covariate
}
\description{
In cape, covariates are coded as genomic loci. To use a phenotype (e.g. male/female) or experimental treatment (e.g. case/control) as a covariate, this factor must be moved to the genotype matrix. \code{create.covar} removes the specified phenotype from the phenotype matrix and places it in the genotype matrix. It is assigned to chromosome 0. 
}
\usage{
create.covar(data.obj, pheno.which)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{pheno.which}{
	Either a column number or a character string indicating which phenotype should be transferred to the genotype matrix. If a character string is used, it must exactly identify the phenotype column header of interest.
}
}
\value{
This function returns the data object with the designated phenotype removed from the phenotype matrix and added to the genotype matrix
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
\examples{
data(obesity.cross)
str(obesity.cross)
obesity.cross <- create.covar(obesity.cross, pheno.which = "mom")
str(obesity.cross)
}
\keyword{manip}