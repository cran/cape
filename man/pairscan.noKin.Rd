\name{pairscan.noKin}
\alias{pairscan.noKin}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Performs a pairscan without a kinship correction.
	}
\description{
This function performs the pairwise regression on all selected marker pairs. The phenotypes used can be either eigentraits or raw phenotypes.
}
\usage{
pairscan.noKin(data.obj, pheno.mat, geno.mat, 
covar.table, paired.markers, n.perm, verbose,
n.cores = NULL)
}
%- maybe also usage for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{pheno.mat}{
	A matrix holding the phenotypes in columns and individuals in rows
	}
	\item{geno.mat}{
	A matrix of genotype values with markers in columns and individuals in rows. Elements indicate the probability of an alternate allele at each marker for each individual.
	}
  \item{covar.table}{
	A matrix holding values for covariates in columns and individuals in rows.
	}
	\item{paired.markers}{
		A two-column matrix indicating which marker pairs to test in the pairscan.
	}
	\item{n.perm}{
		The number of total permutations to be performed.
		}
  \item{verbose}{
	A logical value indicating whether the progress of the scan should be printed to the screen.
	}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing. If NULL, the choice is made automatically.
	}
}
\value{
This function adds an element to the pairscan object reporting the results of the pair-wise scan:
\item{pairscan.effects}{A matrix with five columns indicating the names of markers 1 and 2, their effect sizes, and the effect size of their interaction}
\item{pairscan.se}{A table of the standard errors from the test on each marker pair. The columns are identical to those described for pairscan.effects}
\item{model.covariance}{This is a table in which each row is the linearized matrix of the variance-covariance matrix of each pairwise regression.}
}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{pairscan}}, \code{\link{pairscan.kin}}
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}

