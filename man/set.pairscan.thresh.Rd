\name{set.pairscan.thresh}
\alias{set.pairscan.thresh}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Manually set the threshold by which markers are selected for the pair scan
}
\description{
This function sets the threshold used to filter genetic markers for the pair scan. The threshold is a numeric t statistic (\eqn{\beta/\sigma}). If use of the threshold is indicated in \code{\link{select.markers.for.pairscan}}, only markers with t statistics exceeding the threshold will be selected for the pairscan. Thresholding is useful when the number of genotyped markers is large and an exhaustive pairwise scan is impractical.
}
\usage{
set.pairscan.thresh(data.obj, pairscan.thresh)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{pairscan.thresh}{
		A numeric value indicating the t statistic (\eqn{\beta/\sigma}) to be used as the threshold.
		}
}
\value{
This function returns the data object containing the specified pairscan threshold.
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{select.markers.for.pairscan}}, \code{\link{pairscan}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{manip}

