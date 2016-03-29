\name{sortByThenBy}
\alias{sortByThenBy}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Write out a cape data object to .csv format.
}
\description{
This function sorts a table by two different columns.
}
\usage{
sortByThenBy(tableX, sort.cols = c(1,2), 
col.type = c("c", "n"), decreasing = FALSE, 
return.order = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tableX}{
	The matrix to be sorted.
}
  \item{sort.cols}{
	A two-element numeric vector indicating which columns should be used for sorting, and in which order.
}
	\item{col.type}{
	A two-element character vector indicating whether the two sorting columns contain character values ("c") or numeric values ("n")
	}
	\item{decreasing}{
		A logical value indicating whether the values should be sorted in decreasing order (TRUE) or increasing order (FALSE).
		}
	\item{return.order}{
		If this value is FALSE, the sorted table is returned. If TRUE, a two-column matrix is returned indicating the sorting order for the table. To order the initial table, sort first by the first column, and then by the second. This table can be useful if more than one table is to be sorted in the same order.
	}
}
\value{
This function returns either a sorted table (when return.order is TRUE), or a table indicating the row order of the sorted table (when return.order is FALSE).
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

