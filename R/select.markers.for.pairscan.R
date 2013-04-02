select.markers.for.pairscan <-
function(data.obj, use.pairs.threshold = FALSE, specific.markers = NULL){
	
	
	#take our various parts of the data object
	#for clarity of code later on
	scanone.result <- data.obj$singlescan.results
	covar.flags <- data.obj$covar.flags
	
	#we need to use the covariates determined from the 1D scan
	#in the 2D scan, so stop if the 1D scan has not been performed
	if(length(scanone.result) == 0){
		stop("singlescan() must be run before selecting markers for pairscan\n")
		}	

	
	#if we are using all the possible pairs, take out the whole
	#genotype matrix, otherwise, threshold by the pairscan.threshold
	#we also need to change the covar flags matrix
		if(use.pairs.threshold){
			pairs.thresh <- data.obj$pairscan.thresh
			above.thresh <- unique(sort(unlist(lapply(scanone.result, function(x) which(x[,"t.stat"] >= pairs.thresh)))))
			pair.geno <- data.obj$geno[,above.thresh]
			pair.covar.flags <- covar.flags[above.thresh,]
			}else{
				pair.geno <- data.obj$geno
				pair.covar.flags <- covar.flags
				}

		data.obj$geno.for.pairscan <- pair.geno
		data.obj$covar.for.pairscan <- pair.covar.flags	

		#make sure that the genotype matrix is linearly independent
		geno.ind <- get.linearly.independent(data.obj)
		geno <- geno.ind[[1]]
		rejected.markers <- geno.ind[[2]]
		if(length(rejected.markers) > 0){
			pair.covar.flags <- pair.covar.flags[match(colnames(geno), rownames(pair.covar.flags)),]
			message("\n", length(rejected.markers), " marker(s) rejected due to linear non-independence.\n For more information see markers.removed.txt")
			write.table((data.obj$marker.names)[sort(rejected.markers)], "markers.removed.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
			}
					
	
		#replace the geno and covar flags objects in data.obj with the 
		#linear independent versions
		#sort the genotype matrix and covariate flags before
		#returning them.
		
		marker.order <- match(colnames(geno), colnames(data.obj$geno))
		geno <- geno[,order(marker.order)]
		data.obj$geno.for.pairscan <- geno
		
		marker.order <- match(rownames(pair.covar.flags), colnames(data.obj$geno))
		covar.flags <- pair.covar.flags[order(marker.order),]
		data.obj$covar.for.pairscan <- covar.flags
		# identical(rownames(covar.flags), colnames(geno))

		if(!is.null(specific.markers)){
			markers.which <- get.col.num(geno, specific.markers)
			new.geno <- geno[, markers.which]
			new.covar <- covar.flags[markers.which, ]
			data.obj$geno.for.pairscan <- new.geno
			data.obj$covar.for.pairscan <- new.covar
			}


		return(data.obj)
	
}
