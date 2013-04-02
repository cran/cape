\name{rz.transform}
\alias{rz.transform}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Normalize a vector using rank normalization
}
\description{
This function normalizes numeric vectors to satisfy the assumptions of linear regression. Rank normalization is used to fit the values in the original vector to a normal distribution with a standard deviation of 1.
}
\usage{
rz.transform(x, jitter = FALSE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{x}{
A numeric vector to be normalized
}
  \item{jitter}{
A logical value indicating how rank ties should be handled. If FALSE, ties are not broken. If TRUE, a small jitter factor is added to each term to break ties.
}
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
\seealso{
\code{\link{rank}}
}
\examples{
data(obesity.cross)
phenotype <- obesity.cross$pheno[,"glucose"]
norm.pheno <- rz.transform(phenotype)
hist(phenotype, main = "Histogram of Original\nGlucose 
Distribution", xlab = "raw glucose")
hist(norm.pheno, main = "Histogram of Normalized\nGlucose 
Distribution", xlab = "normalized glucose")
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}
\keyword{distribution}% __ONLY ONE__ keyword per line
