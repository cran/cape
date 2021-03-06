% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_color2.R
\name{get_color2}
\alias{get_color2}
\title{Generate color ramp}
\usage{
get_color2(col_name, col_gap = 10, test_plot = FALSE)
}
\arguments{
\item{col_name}{A string indicating a color for the ramp}

\item{col_gap}{A number indicating the size of the gap
between colors}

\item{test_plot}{Whether to plot colors returned by the function}
}
\value{
A vector of colors forming a color ramp based on the 
input name and numeric gap.
}
\description{
This function generates color ramps based on an input
color. The colors will range from two shades lighter
to two shades darker than the original color, based 
on the color gap specified
}
\keyword{internal}
