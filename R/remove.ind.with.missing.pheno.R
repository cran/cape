remove.ind.with.missing.pheno <-
function(data.obj){
	
	pheno <- data.obj$pheno
	geno <- data.obj$geno
	
	#make sure all the phenotypes are numeric
	pheno <- apply(pheno, 2, as.numeric)
	
	#see if there are any individuals with missing 
	#phenotypes data. Remove individuals with 
	#missing data from both geno and pheno matrices.

	na.rows <- which(is.na(pheno), arr.ind = TRUE)
	rows.to.remove <- sort(unique(na.rows[,1]))

	#stop and warn if all individuals are about to be removed
	if(length(rows.to.remove) == dim(pheno)[1]){
		stop("All individuals have missing phenotype data. Try narrowing your phenotype selection with select.pheno()")
		}
	
	if(length(rows.to.remove) > 0){
		
		#take out the missing entries from the phenotype object
		final.pheno <- pheno[-rows.to.remove,]

		#check the dimensions of the genotype object and remove
		#individuals accordingly

		final.geno <- geno[-rows.to.remove,]

		if(!is.null(data.obj$raw.pheno)){
			data.obj$raw.pheno <- data.obj$raw.pheno[-rows.to.remove,]
			}

		message("\nThere were ", length(rows.to.remove), " individuals removed due to missing phenotypes.")

		}else{

			#if there are no rows to remove, just return the original arrays
			final.pheno <- pheno
			final.geno <- geno

			}
	
	data.obj$pheno <- final.pheno
	data.obj$geno <- final.geno

	return(data.obj)
	
}
