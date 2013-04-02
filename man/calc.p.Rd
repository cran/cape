\name{calc.p}
\alias{calc.p}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Calculate P Values for Interactions Based on Permutations
}
\description{
This function uses the perumtation results to calculate empirical p values for the variant-to-variant influences calculated by \code{\link{error.prop}}. It also adjusts these p values using Holm's step-down procedure.
}
\usage{
calc.p(data.obj)
}

\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
}
\value{
The data object is returned with a new list with two elements. The elements correspond to the two directions of influence: marker1 to marker2 and marker2 to marker1. Each element contains a table with the source and target variants, the empirical p values, and the adjusted p values, along with the effect size, standard error and t statistic for each interaction.
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
\examples{
# not run
# obesity.cross <- calc.p(obesity.cross)
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{models}% __ONLY ONE__ keyword per line
