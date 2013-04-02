get.eigentraits <-
function(data.obj, scale.pheno = TRUE, normalize.pheno = TRUE){

	#first make sure there are no individuals
	#with missing phenotypes. This also makes 
	#sure the phenotypes are numeric
	data.obj <- remove.ind.with.missing.pheno(data.obj)

	pheno <- data.obj$pheno

	#if there are fewer phenotypes than individuals...
	if(dim(pheno)[1] < dim(pheno)[2]){
			message("\nThere must be more individuals than phenotypes.\nPlease select fewer phenotypes to analyze or\nperform a dimension reduction on the\nphenotypes before reading in the data.")
			return(data.obj)
			}

	#This function mean centers and standardizes a vector
	center.std <- function(v){
		mean.v <- mean(v)
		centered <- v - mean.v
		sd.v <- sd(v)
		final.v <- centered/sd.v
		return(final.v)
		}

	if(normalize.pheno){
		pheno <- apply(pheno, 2, rz.transform)
		}

	if(scale.pheno){
		pheno <- apply(pheno, 2, center.std) #mean center and standardize the phenotypes
		}


	data.obj$pheno <- pheno #replace the raw phenotypes with scaled, normalized phenotypes (if we have done those things)

	svd.pheno <- svd(pheno)
	

	#add the eigentraits and singular values to the data object
	data.obj$ET <- svd.pheno$u
	data.obj$singular.values <- svd.pheno$d
	data.obj$right.singular.vectors <- svd.pheno$v
	colnames(data.obj$ET) <- paste("ET", 1:length(data.obj$pheno[1,]), sep = "")

	
	return(data.obj)
	
}
