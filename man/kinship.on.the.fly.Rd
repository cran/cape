\name{kinship.on.the.fly}
\alias{kinship.on.the.fly}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Kinship correction for genotype and phenotype.
}
\description{
This function calculates corrected genotype, phenotype, and covariate values using previously calculated kinship matrices.
}
\usage{
kinship.on.the.fly(kin.mats, geno, chr1 = NULL, 
chr2 = NULL, phenoV, covarV = NULL)
}

\arguments{
  \item{kin.mats}{
	The list of kinship matrices calculated by \code{\link{kinship}}.
	}
\item{geno}{
	The genotype matrix.
	}
\item{chr1}{
	Kinship matrices are calculated in a leave-one-out or leave-two-out manner, i.e. the chromosome(s) on which the markers being tested are located are left out of the kinship correction. This argument indicates the chromosome on which one marker of the marker pair resides. It indicates which kinship matrix to use for the data correction. If chr2 == chr1, the leave-one-out kinship matrix for that chromosome will be used. Otherwise, the leave-two-out kinship matrix is used. If both chr1 and chr2 are NULL, the kinship matrix for the entire genome will be used.

	}
\item{chr2}{
	Kinship matrices are calculated in a leave-one-out or leave-two-out manner, i.e. the chromosome(s) on which the markers being tested are located are left out of the kinship correction. This argument indicates the chromosome on which one marker of the marker pair resides. It indicates which kinship matrix to use for the data correction. If chr2 == chr1, the leave-one-out kinship matrix for that chromosome will be used. Otherwise, the leave-two-out kinship matrix is used. If both chr1 and chr2 are NULL, the kinship matrix for the entire genome will be used.
	}
\item{phenoV}{
	The phenotype vector to be corrected by the appropriate kinship matrix.
	}
\item{covarV}{
	The covariate vector to be corrected by the appropriate kinship matrix.
	}
}
\value{
Returns a list of data corrected using the appropriate kinship matrix.
The elements of the list are as follows:
\item{err.cov}{The error covariance matrix from the calculation.}
\item{corrected.pheno}{The kinship-corrected phenotype values.}
\item{corrected.geno}{The kinship-corrected genotype values.}
\item{corrected.covar}{The kinship-corrected covariate values.}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}
\seealso{
\code{\link{kinship}}
}