\name{get.chr.markers}
\alias{get.chr.markers}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Calculate a genome-wide significance threshold for the single-variant scan
}
\description{
This function finds all marker pairs on a given pair of chromosomes
}
\usage{
get.chr.markers(chr.pair, marker.chr, paired.markers)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{chr.pair}{
	The chromosome pair in question.
	}
\item{marker.chr}{
	A matrix with two columns with an entry for each marker pair indicating which chromosome each marker is on.
	}
  \item{paired.markers}{
	A matrix with two columns indicating all marker pairs being tested.
	}
}
\value{
Returns all marker pairs on the chromosomes specified.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

