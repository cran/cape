\name{norm.pheno}
\alias{norm.pheno}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Normalize and mean center phenotypes
}
\description{
This function performs quantile normalization on phenotypes and optionally mean centers them.
}
\usage{
norm.pheno(data.obj, mean.center = TRUE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{mean.center}{
		A logical value. If TRUE the phenotypes are mean centered in addition to being normalized. if FALSE, the phenotypes are not mean centered.
		}
}
\details{
In quantile normalization the values of the phenotype are sorted and replaced with a corresponding value drawn from a normal distribution with the same standard deviation and mean as the original distribution.
Mean centering subtracts the mean phenotype value from each phenotype value yielding a distribution centered around 0.
}
\value{
This function returns data.obj with the normalized phenotypes in place of the original phenotypes. 
}
\note{
 Both normalization and mean centering are highly recommended before obtaining the eigentraits with singular value decomposition (SVD) (see \code{\link{get.eigentraits}}). 
 This normalization procedure can also be performed by \code{\link{get.eigentraits}} 
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{get.eigentraits}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{arith}
\keyword{manip}

