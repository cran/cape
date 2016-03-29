\name{direct.influence}
\alias{direct.influence}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate the significance of direct influences of variant pairs on phenotypes
}
\description{
This function recasts the variant-to-eigentrait effects in terms of variant-to-phenotype effects. It multiplies the \eqn{\beta}-coefficient matrices of each variant (i) and each phenotype (j) (\eqn{\beta_{i}^{j}}) by the singular value matrices (\eqn{V \cdot W^{T}}) obtained from the singular value decomposition performed in \code{\link{get.eigentraits}}. \eqn{\beta_{i}^{j} = V \cdot W^{T}}. It also uses the permutation data from the pairwise scan (\code{\link{pairscan}}) to calculate an empirical p value for the influence of each marker pair on each phenotype. The empirical p values are then adjusted for multiple testing using Holm's step-down procedure.
}
\usage{
direct.influence(data.obj, pairscan.obj, transform.to.phenospace = TRUE, 
pval.correction = c("holm", "fdr", "lfdr"), verbose = FALSE, 
save.permutations = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.obj}{
	The object in which all results are stored. See \code{\link{read.population}}.
	}
	\item{pairscan.obj}{
	The object in which the results from the pairscan are stored. See \code{\link{pairscan}}.
	}
	\item{transform.to.phenospace}{
		A logical value. If TRUE, the influence of each marker on each eigentrait is transformed to the influence of each marker on each of the original phenotypes. If FALSE, no
		transformation is made. If the pair scan was done on eigentraits, the influence of each marker on each eigentrait is calculated. If the pair scan was done on raw phenotypes, 
		the influence of each marker on each phenotype is calculated. The default behavior is to transform variant influences on eigentraits to variant influences on phenotypes.
	}
	\item{pval.correction}{
	One of "holm", "fdr", or "lfdr", indicating whether the p value correction method used should be the Holm step-down procedure, false discovery rate or local false discovery rate respectively.
	}
	\item{verbose}{
		A logical value. If TRUE, the progress of the function is printed to the screen. If FALSE, the default, nothing is printed.
		}
	\item{save.permutations}{
		A logical value indicating whether the data from permutations should be saved. Saving the permutations requires more memory but can be helpful in diagnostics. If save.permutations is TRUE all permutation data are saved in an object called "permutation.data.RData".
	}
}

\value{
This function adds a total of six objects to the data object.
First, a flag (\code{transform.to.phenospace}) is added to the data object to indicate whether variant influences were transformed to phenotype space.

The results from the pairwise scan and the permutations of the pairwise scan are converted to phenospace if specified. These actions each add one object each to the data object (\code{var.to.pheno.influence} and 
\code{var.to.pheno.influence.perm}). Each element is itself a list of matrices corresponding to the original phenotypes. Each matrix contains one row per marker pair (or permutation of a marker pair) and contains the influence coefficient and standard error of the influence coefficient for each pair.

After the coefficients have been transformed to phenotype space, each marker is considered individually and its influence on each phenotype across all marker pair contexts is tabulated. This is done for both the pairwise scan and the permutations of the pairwise scans and adds two new objects (\code{var.to.pheno.test.stat} and \code{var.to.pheno.test.stat.perm}) to the data object. Each object is a list containing one element for each of the original phenotypes. Each element contains a table in which all instances of each marker are listed along with that marker's direct phenotypic influence, the standard error of the influence, and the t statistic (\eqn{\beta/\sigma}) of the influence.

Finally, because each marker can only have one influence on each phenotype, the influences from each marker pair context are filtered to report only the maximum influence of each marker on each phenotype across all marker pair contexts. This process adds an object to the data object called \code{max.var.to.pheno.influence}. This object is a list containing one element per phenotype. It tabulates the maximum influence of each marker on each phenotype, as well as the empirical and Holm's corrected p values associated with each influence.

}
\references{
Carter, G. W., Hays, M., Sherman, A., & Galitski, T. (2012). Use of pleiotropy to model genetic interactions in a population. PLoS genetics, 8(10), e1003010. doi:10.1371/journal.pgen.1003010
Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian journal of statistics, pages 65-70.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{arith}

