\name{select.eigentraits}
\alias{select.eigentraits}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Select a subset of the eigentraits for further analysis
}
\description{
This function selects the specified eigentraits for further analysis. After the singular value decomposition (SVD) of multiple traits by \code{\link{get.eigentraits}}, a subset of the eigentraits, for example the first two, may carry useful information, while others may be dominated by noise. In this case, this function can be used to select the first two eigentraits for use in the analysis. 
}
\usage{
select.eigentraits(data.obj, traits.which = c(1, 2))
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{traits.which}{
A numeric vector indicating which eigentraits should be retained for further analysis.
}
}
\details{
Before the selection of eigentraits, the eigentraits should be examined using \code{\link{plotSVD}}
}
\value{
This function returns the data object with only the selected eigentraits.
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{get.eigentraits}}, \code{\link{plotSVD}}
}
\examples{
# not run
# obesity.cross <- select.eigentraits(obesity.cross, c(1,2))
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}
