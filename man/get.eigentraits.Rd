\name{get.eigentraits}
\alias{get.eigentraits}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Calculate eigentraits from phenotype matrix
}
\description{
This function performs the singular value decomposition (SVD) on the phenotype matrix after first removing individuals with missing data.
The eigentraits are the left singular vectors of the decomposition. This function optionally mean centers and normalizes the phenotype matrix before performing the SVD.
}
\usage{
get.eigentraits(data.obj, scale.pheno = TRUE, normalize.pheno = TRUE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{scale.pheno}{
	A logical value specifying whether the phenotypes should be mean centered before the SVD is performed. The default, and recommended, value is TRUE.
}
  \item{normalize.pheno}{
	A logical value specifying whether the phenotypes should be quantile normalized before the SVD is performed.
}
}
\value{
This function adds three new elements to the data.obj list.
\item{ET}{The left singular vectors from the SVD. These are the eigentraits.}
\item{singular.values}{The singular values from the SVD. These are used later internally to convert variant effects from eigentrait space to phenotype space.}
\item{right.singular.vectors}{The right singular vectors from the SVD. These are used later internally to convert variant effects from eigentrait space to phenotype space.}
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
\note{
There must be more individuals than phenotypes to perform this calculation. An error results if there are more phenotypes than individuals.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{norm.pheno}}, \code{\link{plotSVD}}
}
\examples{
# not run
# obesity.cross <- get.eigentraits(obesity.cross)
}

% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{algebra}

