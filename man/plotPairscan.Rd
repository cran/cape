\name{plotPairscan}
\alias{plotPairscan}
%- Also NEED an \alias for EACH other topic documented here.
\title{ plot the results from \code{\link{pairscan}}
}
\description{
This function plots the results from the pairwise regressions. The matrices are plotted with results coded on a yellow to blue scale. Results from all phenotypes are plotted on the same scale.
}
\usage{
plotPairscan(data.obj, phenotype = NULL, standardized = FALSE, 
show.marker.labels = FALSE, show.chr = TRUE, 
label.chr = TRUE, pdf.label = "Pairscan.Regression.pdf", 
verbose = FALSE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
	\item{phenotype}{
		A character vector indicating which phenotypes should be plotted. All phenotypes are plotted if phenotype is set to NULL.
	}
	\item{standardized}{
		A logical value. If FALSE, the plotted values are the beta coefficients of the interaction in the pairwise model. If TRUE, the t statistics (\eqn{\beta/\sigma}) are plotted.
	 }
	 \item{show.marker.labels}{
	 	A logical value indicating whether marker labels should appear on the axes.
	 }
	 \item{show.chr}{
	 	A logical value. If TRUE, chromosomes are indicated by alternating gray and white blocks along plot axes.
	 	}
	 \item{label.chr}{
	 	A logical value. If TRUE and if show.chr is TRUE, chromosome numbers are printed inside each gray and white block.
	 	}
	 \item{pdf.label}{
	 	A character string indicating the name of the file that the plots should be printed to. When multiple phenotypes are plotted, all penotypes are plotted on the same scale and each is plotted in a different page of the pdf.
	 }
	 \item{verbose}{
	 	A logical value indicating whether the progress of the plotting function should be reported. Large data sets can take a long time to process and plot.
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
