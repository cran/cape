\name{bin.vector}
\alias{bin.vector}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Bin a continuous vector into discrete values.
}
\description{
This function bins a vector with continuous values into discrete values by assigning each continuous value to its neares discrete value.
}
\usage{
bin.vector(vectorX, bins = seq(0,1,0.5))
}

\arguments{
  \item{vectorX}{
	The vector to be binned
	}
\item{bins}{
	The discrete values into which vectorX should be binned.
	}
}
\value{
Returns a binned vector with only the discrete values specified by bins.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}% __ONLY ONE__ keyword per line
