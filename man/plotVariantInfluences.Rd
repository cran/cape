\name{plotVariantInfluences}
\alias{plotVariantInfluences}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot variant-to-variant influences
}
\description{
This function plots the reparameterized influences of variants on each other. The epistatic interactions from the pairwise scan are reparameteriezed to the terms \eqn{m_{12}} and \eqn{m_{21}}, where the subscripts indicate the source and target variants respectively. These terms are interpreted as the effect that the source variant exerts on the target variant when both are present. Negative influences represent suppression while positive influences represent enhancement.
}
\usage{
plotVariantInfluences(data.obj, p.or.q = 0.05, 
all.markers = FALSE, standardize = FALSE, 
not.tested.col = "lightgray", show.marker.labels = FALSE, 
show.chr = TRUE, label.chr = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
	\item{p.or.q}{
		A threshold indicating the maximum adjusted p value considered significant. If an fdr method has been used to correct for multiple testing, this value specifies the maximum q value considered significant. Only marker pairs with p or q values below this threshold will be plotted.
	}
	\item{all.markers}{
		A logical value. If TRUE all markers are plotted. If FALSE only markers tested in the pair scan are plotted.
}
	\item{standardize}{
		A logical value. If TRUE the values in each entry of the plotted matrix are the standardized effect sizes \eqn{\beta/\sigma}. If FALSE, the raw \eqn{\beta} values are plotted.
	}
	\item{not.tested.col}{
		A color name used to mark variant pairs that were not tested due to linkage. The color distinguishes these pairs from those that were tested, but do not have significant interactions. The default color is light "gray". This can be changed to "white" or FALSE if no marking is desired.
	}
	\item{show.marker.labels}{
		A logical value indicating whether the marker labels should be printed along the plot axes.
	}
	 \item{show.chr}{
	 	A logical value. If TRUE, chromosomes are indicated by alternating gray and white blocks along plot axes.
	 	}
	 \item{label.chr}{
	 	A logical value. If TRUE and if show.chr is TRUE, chromosome numbers are printed inside each gray and white block.
	 	}
}
\value{
No values are returned. The function generates a plot of the reparameterized source to target effects.
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}

\seealso{
\code{\link{pairscan}}}
\examples{
# not run
# plotVariantInfluences(obesity.cross)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{hplot}

