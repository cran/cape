\name{plotPairscan}
\alias{plotPairscan}
%- Also NEED an \alias for EACH other topic documented here.
\title{ plot the results from \code{\link{pairscan}}
}
\description{
This function plots the results from the pairwise regressions. The matrices are plotted with results coded on a yellow to blue scale. Results from all phenotypes are plotted on the same scale.
}
\usage{
plotPairscan(data.obj, phenotype = "ET1", standardized = FALSE, 
pdf.label = "Pairscan.Regression.pdf")
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
	\item{phenotype}{
		A character vector indicating which phenotypes should be plotted. All phenotypes can be plotted by setting phenotype to NULL.
	}
	\item{standardized}{
		A logical value. If FALSE, the plotted values are the beta coefficients of the interaction in the pairwise model. If TRUE, the t statistics (\eqn{\beta/\sigma}) are plotted.
	 }
	 \item{pdf.label}{
	 	A character string indicating the name of the file that the plots should be printed to. When multiple phenotypes are plotted, all penotypes are plotted on the same scale and each is plotted in a different page of the pdf.
	 }
}
\value{
No values are returned.
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{pairscan}}
}
\examples{
# not run
# plotPairscan(data.obj)
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{hplot}
