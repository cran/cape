\name{plotSVD}
\alias{plotSVD}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Plot the results of the singular value decomposition of the phenotype matrix
}
\description{
This function plots the results of the singular value decomposition (SVD) of the phenotype matrix. This plot is useful for selecting eigentraits for further analysis.
}
\usage{
plotSVD(data.obj, orientation = c("vertical", "horizontal"), 
neg.col = "blue", pos.col =  "brown", light.dark = "dark")
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{orientation}{
	A character string ("vertical", or "horizontal") indicating whether the plot should be vertically or horizontally oriented.
}
	 \item{neg.col}{
	 	The color in which to plot the negative values. Colors can be one of the following: "green", "purple", "orange", "blue", "brown", "yellow", or "gray."
	 }
	 \item{pos.col}{
	 	The color in which to plot the positive values. Colors can be one of the following: "green", "purple", "orange", "blue", "brown", "yellow", or "gray."
	 }
	 \item{light.dark}{
	 	Either "light" or "dark" depending on the desired range of colors used.
	 }
}
\value{
Invisibly returns the fraction of variance accounted for by each ET.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{get.eigentraits}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{hplot}

