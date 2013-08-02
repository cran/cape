\name{segment.region}
\alias{segment.region}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Divide a region into equal parts.
}
\description{
This function is useful for placing a given number of elements into a constrained space on a plot. It returns coordinates for the desired number of points in the given region. The points can either be centered within the region or flush to the ends.
}
\usage{
segment.region(region.min, region.max, num.points, alignment = c("center", "ends"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{region.min}{
	A numerical value indicating the minimum value of the region.
}
  \item{region.max}{
	A numerical value indicating the maximum value of the region.
}
  \item{num.points}{
	The number of points to place in the region.
}
  \item{alignment}{
	A character element indicating whether the points should be centered within the region or whether they should extend to the ends of the region.
}
}
\value{
A vector of coordinates is returned specifying the coordinates at which points should be plotted to be evenly distributed in the region.
}


\examples{
num.points <- 5
coords.1 <- segment.region(region.min = 5, region.max = 7, 
num.points = num.points, alignment = c("ends"))
coords.2 <- segment.region(region.min = 5, region.max = 7, 
num.points = num.points, alignment = c("center"))
plot(x = coords.1, y = rep(0.25, num.points), xlim = c(4,8), 
ylim = c(0,1),col = "blue", pch = 16, xlab = "", ylab = "", 
axes = FALSE)
points(x = coords.2, y = rep(0.75, num.points), xlim = c(4,8), 
ylim = c(0,1),col = "purple", pch = 16)
abline(v = 5); abline(v = 7); axis(1)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{internal}

