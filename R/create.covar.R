create.covar <-
function(data.obj, pheno.which){

	pheno <- data.obj$pheno
	geno <- data.obj$geno
	chromosome <- data.obj$chromosome
	marker.location <- data.obj$marker.location
	
	new.marker.locale <- get.col.num(pheno, pheno.which)
	
	if(length(new.marker.locale) == 0){
		return(data.obj)
		}
	
	if(length(unique(new.marker.locale)) < length(new.marker.locale)){
		stop("Phenotype labels cannot be substrings of other phenotypes.")
		}
	
	new.marker <- matrix(pheno[,new.marker.locale], ncol = length(new.marker.locale))

	new.geno <- cbind(geno, new.marker)
	
	
	colnames(new.geno) <- 1:dim(new.geno)[2]
	rownames(new.geno) <- rownames(geno)
		
	new.marker.names <- c(data.obj$marker.names, pheno.which)
	chromosome <- c(chromosome, rep(0, length(pheno.which)))
	marker.location <- c(marker.location, 1:length(pheno.which))
	
	data.obj$geno <- new.geno
	data.obj$chromosome <- chromosome
	data.obj$marker.location <- marker.location
	data.obj$marker.names <- new.marker.names
	
	#take the phenotypes made into markers out of the phenotype matrix
	new.pheno <- pheno[,-new.marker.locale]
	data.obj$pheno <- new.pheno

	return(data.obj)
	
	}
