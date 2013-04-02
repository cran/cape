\name{select.by.ind}
\alias{select.by.ind}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
Subset a cross object to include specific individuals
}
\description{
This function subsets a cross to include individuals based on either phenotypic or genotypic values. For example, this function can subset a cross to include all individuals with a phenotype value greater than x or a genotype value equal to y.
}
\usage{
select.by.ind(data.obj, geno.or.pheno = pheno, expr)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{geno.or.pheno}{
A character value, either "geno", or "pheno" to specify which matrix should be used to subset individuals. 
}
  \item{expr}{
An quoted expression used to subset individuals.
}
}
\value{
The cross object is returned including only the individuals meeting the criteria in the provided expression.
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{select.by.chr}}
}
\examples{
data(obesity.cross)
hist(obesity.cross$pheno[,"insulin"], main = "original insulin distribution",
xlab = "insulin (ng/mL)", xlim = c(0, 30))
obesity.cross <- select.by.ind(obesity.cross, "pheno", "insulin < 25")
hist(obesity.cross$pheno[,"insulin"], main = "subset insulin distribution", 
xlab = "insulin (ng/mL)", xlim = c(0, 30))
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}
