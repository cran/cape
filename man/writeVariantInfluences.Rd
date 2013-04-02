\name{writeVariantInfluences}
\alias{writeVariantInfluences}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write the final results to a file
}
\description{
This function takes in the final data object and writes the variant influences that are at or below the specified significance level
to a file in the current working directory.
}
\usage{
writeVariantInfluences(data.obj, p.or.q = 0.05, 
filename = "Variant.Influences.csv", delim = ",")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
}
  \item{p.or.q}{
	A threshold indicating the maximum adjusted p value considered significant. If an fdr method has been used to correct for multiple testing, this value specifies the maximum q value considered significant. Only marker pairs with p or q values below this threshold will be plotted.
}
  \item{filename}{
	A character vector specifying the name of the file.
	}
	\item{delim}{
		A character string indicating the delimeter in the data file. The default indicates a comma-separated file (",").
	}
}
\value{
This function writes a table of direct influences to a file. It also returns the table invisibly, i.e. if the output of the function is assigned to an object, the object will contain the table of influences. Otherwise, nothing is returned.
}
\author{
Anna L. Tyler, Wei Lu, Justin J. Hendrick, Vivek M. Philip, and Gregory W. Carter
}
\examples{
# not run
# here the table is written to a file, but nothing is returned.
# writeVariantInfluences(obesity.cross)

# here the table is written to a file, and returned to the 
# object sig.table
# sig.table <- writeVariantInfluences(obesity.cross)
# print(sig.table)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{IO}

