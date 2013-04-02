select.by.chr <-
function(data.obj, chr, include.covariates = TRUE){
	
	if(include.covariates){
		chr <- c(chr, 0)
		}

	geno <- data.obj$geno
	
	chr.list <- data.obj$chromosome
	chr.locale <- which(chr.list %in% chr)
	
	sub.geno <- geno[,chr.locale]
	data.obj$geno <- sub.geno
	data.obj$chromosome <- data.obj$chromosome[chr.locale]
	data.obj$marker.location <- data.obj$marker.location[chr.locale]

	return(data.obj)
	
	}
