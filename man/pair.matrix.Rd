\name{pair.matrix}
\alias{pair.matrix}
%- Also NEED an \alias for EACH other topic documented here.
\title{ Given a vector of elements, create a two-column matrix listing all pairs of elements
}
\description{
This function takes in a vector and returns a matrix of all pairs of the elements. The order of pairs can be ignored or taken into account, and self-pairs can either be included or ignored.
}
\usage{
pair.matrix(elements, ordered = FALSE, self.pairs = FALSE)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{elements}{
	A vector of elements for which all pairs are wanted. This vector can include numeric values, strings, or a mixture of both.
}
  \item{ordered}{
	A logical value indicating whether order of pairs is important. TRUE indicates that pair 'a-b' is different from 'b-a'. A value of FALSE will return only one of these pairs.
}
  \item{self.pairs}{
	A logical value indicating whether elements should be paired with themselves ('a-a'). A value of TRUE will return a matrix including self pairs.
}
}
\value{
	A two-column matrix of all pairwise combinations of the vector elements is returned.
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\examples{
el <- letters[1:4]
letter.pairs <- pair.matrix(el, ordered = FALSE, self.pairs = FALSE)
print(letter.pairs)

el <- 1:4
num.pairs <- pair.matrix(el, ordered = FALSE, self.pairs = TRUE)
print(num.pairs)
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

