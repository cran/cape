\name{calc.p}
\alias{calc.p}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Calculate P Values for Interactions Based on Permutations
}
\description{
This function uses the perumtation results to calculate empirical p values for the variant-to-variant influences calculated by \code{\link{error.prop}}. It can also optionally adjust these p values using Holm's step-down procedure, false discovery rate (fdr), or local false discovery rate (lfdr).
}
\usage{
calc.p(data.obj, pval.correction = c("holm", "fdr", "lfdr", "none"))
}

\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{pval.correction}{
		One of "holm", "fdr", "lfdr" or "none", indicating whether the p value correction method used should be the Holm step-down procedure, false discovery rate, local false discovery, or no correction rate respectively.
	}
}
\value{
The data object is returned with a new list with two elements. The elements correspond to the two directions of influence: marker1 to marker2 and marker2 to marker1. Each element contains a table with the source and target variants, the empirical p values, and the adjusted p values, along with the effect size, standard error and t statistic for each interaction.
}
\references{
Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian journal of statistics, pages 65-70.
Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Societ. Series B (Methodological), 289-300.
Liao, J.G., Lin, Y., Selvanayagam, Z.E., & Shih, W.J. (2004). A mixture model for estimating the local false discovery rate in DNA microarray analysis. Bioinformatics, 20(16), 2694-2701. doi:10.1093/bioinformatics/bth310
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{models}% __ONLY ONE__ keyword per line
