\name{kinship}
\alias{kinship}
%- Also NEED an \alias for EACH other topic documented here.
\title{ 
	Calculate all leave-one-out or leave-two-out kinship matrices.
}
\description{
This function calculates kinship matrices for all individual or pairs of chromosomes being tested.
}
\usage{
kinship(geno, all.chr, chr.which = NULL, 
chr.pairs = FALSE, sample.kinship = FALSE, 
num.samples = 100, n.per.sample = 10, 
run.parallel = TRUE, add.full.kin = FALSE, 
verbose = FALSE, n.cores = NULL)
}

\arguments{
  \item{geno}{
	The genotype matrix.
	}
\item{all.chr}{
	A vector indicating the chromosome on which each marker in the genotype matrix resides.
	}
\item{chr.which}{
	A vector indicating for which chromosomes kinship matrices should be calculated.
	}
\item{chr.pairs}{
	A logical value. If FALSE, kinship corrections are calculated on a leave-one-out basis only. There will be one kinship matrix for each chromosome left out of the calculation. If TRUE, kinship corrections are calculated based on leaving two chromosomes out, and a kinship matrix will be calculated for each pair of chromosomes.
	}
\item{sample.kinship}{
	A logical value indicating whether kinship matrices should be sampled (TRUE) or calculated directly (FALSE). When there are very large numbers of markers in the genotype matrix sampling can approximate the kinship matrix in less time than an exact calculation.
	}
\item{num.samples}{
	If sample.kinship is TRUE, this number indicates how many samples will be used to calculate the kinship matrix for each leave-one-out or leave-two-out instance.
	}
\item{n.per.sample}{
	If sample.kinship is TRUE, this number indicates the number of markers that should be used in each sample of the kinship matrix.
	}
\item{run.parallel}{
	A logical value indicating whether the process should be run in parallel.
	}
\item{add.full.kin}{
	A logical value. If TRUE, a kinship matrix using the entire genome will be calculated. This is used when covariates are being tested. They are corrected based on a kinship matrix calculated from the entire set of markers.
	}
\item{verbose}{
	A logical value. If TRUE, the progress of the process is printed to the screen.
	}
	\item{n.cores}{
		An integer specifying the number of cores to be used in parallel processing. If NULL, the choice is made automatically.
	}
}
\value{
Returns a list of kinship matrices, one for each chromosome or each pair of chromosomes depending on whether chr.pairs is TRUE or FALSE.
}
% Add one or more standard keywords, see file KEYWORDS in the
% R documentation directory.
\keyword{internal}
\seealso{
\code{\link{kinship.on.the.fly}}
}
\references{
  Cheng R, Abney M, Palmer AA and Skol AD (2011) QTLRel: an R Package for Genome-wide Association Studies in which Relatedness is a Concern.
  BMC Genet. 12: 66

  Cheng R, Lim JE, Samocha KE, Sokoloff G, Abney M, Skol AD and Palmer AA (2010) Genome-wide association studies and the problem of
  relatedness among advanced intercross lines and other highly recombinant populations. Genetics 185: 1033-1044

Kang, H.M., Zaitlen, N.A., Wade, C.M., Kirby, A., Heckerman, D., Daly, M. J., and Eskin, E. (2008). Efficitne control of population structure in model organism association mapping. Genetics, 178(3), 1709-1723.
}
