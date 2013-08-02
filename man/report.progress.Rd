\name{report.progress}
\alias{report.progress}
%- Also NEED an \alias for EACH other topic documented here.
\title{
Print the progress of a function to the screen
}
\description{
This is an internal function that takes in the current position in a process and the total length of the process and prints to the screen the percent progress at the specified increments.
}
\usage{
report.progress(current, total, percent.text = 10, percent.dot = 2)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{current}{
	A numeric value specifying the current position in the process.
	}
	\item{total}{
		A numeric value specifying the total number of steps in the process
	}
	\item{percent.text}{
		A numeric value indicating the frequency at which the percent done should be reported. The default value is 10 indicating that text indicating the progress is printed to the screen when the process is 10, 20, 30 percent done, etc.
	}
	\item{percent.dot}{
		A numeric value indicating the frequency at which dots in between text reports are printed. The default value is 2 indicating that a dot is printed to the screen when the process is 2, 4, 6 percent done, etc.
	}
}
\value{
This function does not return any values but prints values to the screen.
}

\examples{
	
#report progress every 10 percent, print dots every 2%
all.steps <- 1:100
for(i in all.steps){
	report.progress(i, total = length(all.steps), percent.text = 10, percent.dot = 2)
}

#report progress every 25 percent, print dots every 5%
all.steps <- 1:100
for(i in all.steps){
	report.progress(i, total = length(all.steps), percent.text = 25, percent.dot = 5)
}

}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{arith}
\keyword{internal}

