\name{obesity.cross}
\alias{obesity.cross}
\docType{data}
\title{ 
Mouse cross data from Reifsnyder et al. (2000)
}
\description{
Data from a cross between non-obese, non-diabetic (NON) mice, and diabetes-prone New Zealand obese (NZO) mice. The experiment is described in Reifsnyder et al. (2000). The data object is a list. The first element is a matrix of phenotype data. In includes insulin levels (ng/mL), plasma glucose levels (mg/dL), total body weight (g), all measured at age 24 weeks. The second element is the genotype matrix containing the genotype of each mouse at each of 85 markers across the  genome. The final elements are the chromosome vector, which indicates which chromosome on which each marker is found, and a vector of chromosomal positions indicating the location of each marker on the chromosome in centimorgans (cM).
}
\usage{data(obesity.cross)}
\format{
  The format is:
List of 4
 $ pheno          : num [1:204, 1:3] 63.3 31 43.3 33.3 35.3 15.9 23.6 NA 16 22 ...
  ..- attr(*, dimnames)=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:3] body_weight sex
 $ geno           : num [1:204, 1:85] 0.5 0.5 0.5 0 0.5 0 0.5 NA 0 0 ...
  ..- attr(*, dimnames)=List of 2
  .. ..$ : NULL
  .. ..$ : chr [1:85] D1Mit296 D1Mit211 D1Mit411 D1Mit123 ...
 $ chromosome     : chr [1:85] 1 1 1 1 ...
 $ marker.location: num [1:85] 2.08 10.59 12.62 17.67 22.88 ...
}
\source{
Reifsnyder, P. C. (2000). Maternal Environment and Genotype Interact to Establish Diabesity in Mice. Genome Research, 10(10), 1568-1578.
}
\keyword{datasets}
