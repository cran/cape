% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bin_vector.R
\name{bin_vector}
\alias{bin_vector}
\title{Snap continuous values to a grid}
\usage{
bin_vector(vectorX, bins = seq(0, 1, 0.5))
}
\arguments{
\item{vectorX}{A vector of numeric values to bin}

\item{bins}{A vector of values to snap the values in vectorX to.}
}
\value{
A vector the same length as vectorX in which 
each value in vectorX has been sent to the nearest value
in bins. For example, if bins is c(0, 0.5, 1), and vectorX
contains a 0.49. That 0.49 value will be sent to 0.5.
The returned vector contains only the values specified in bins.
}
\description{
This function bins a continuously valued vector
based on user-defined bins.
It is useful for binning continuously valued genotypes.
Each value in the matrix gets shifted to the nearest
provided in the argument bins.
}
\keyword{internal}
