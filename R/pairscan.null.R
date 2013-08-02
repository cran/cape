pairscan.null <-
function(data.obj, scan.what = c("eigentraits", "raw.traits"), total.perm = NULL, n.top.markers = 50, verbose = FALSE){

	if(is.null(total.perm)){
		stop("The total number of permutations must be specified.")
		}

	scanone.result <- data.obj$singlescan.results
	
	#we need to use the covariates determined from the 1D scan
	#in the 2D scan, so stop if the 1D scan has not been performed
	if(length(scanone.result) == 0){
		stop("singlescan() must be run on this object before pairscan()\n")
		}	


	#If the user does not specify a scan.what, 
	#default to eigentraits, basically, if eigen,
	#et, or ET are anywhere in the string, use the
	#eigentraits, otherwise, use raw phenotypes
	type.choice <- c(grep("eig", scan.what), grep("ET", scan.what), grep("et", scan.what)) #look for any version of eigen or eigentrait, the user might use.
	if(length(type.choice) > 0){ #if we find any, use the eigentrait matrix
		pheno <- data.obj$ET
		}else{
			pheno <- data.obj$pheno #otherwise, use the raw phenotype matrix
			}


	geno <- data.obj$geno.for.pairscan
	covar.flags <- data.obj$covar.for.pairscan

	if(is.null(geno)){
		stop("select.markers.for.pairscan() must be run before pairscan()")
		}
		

	#make a list to hold the results. 
	#Tables from each of the phenotypes will be
	#put into this list
	results.perm.list <- vector(mode = "list", length = dim(pheno)[2])

	 	 
	#generate the null distribution for the pairscan using the 
	#total number of permutations and the permuted results from
	#singlescan.
	
		for(p in 1:length(scanone.result)){ 

			if(verbose){
				cat("\nGenerating null distribution for ", colnames(pheno)[p], ":\n", sep = "")
				}
			
			final.perm <- 1
			j <- 1
			while(final.perm < total.perm){
				if(verbose){cat(round((final.perm/total.perm)*100), "%...", sep = "")} 
				# if(verbose){report.progress(current = final.perm, total = total.perm)}
				perm.pheno <- pheno[sample(1:dim(pheno)[1], dim(pheno)[1], replace = FALSE),p]
				single.scan.result <- one.singlescan(phenotype.vector = perm.pheno, genotype.mat = geno, covar.vector = covar.flags[,p])
				
				geno.order <- order(single.scan.result[,"t.stat"], decreasing = TRUE)
				top.geno.mat <- geno[,geno.order[1:n.top.markers]]
				top.geno.covar.flags <- covar.flags[geno.order[1:n.top.markers],p]
				top.marker.pairs <- pair.matrix(colnames(top.geno.mat))
				
				pairscan.results <- one.pairscan(phenotype.vector = perm.pheno, genotype.matrix = top.geno.mat, covar.vector = top.geno.covar.flags, pairs.matrix = top.marker.pairs, n.perm = 0, verbose = FALSE)
				
				#integrate the results into the permutation object
				one.perm <- pairscan.results[[1]]
				if(final.perm == 1){ #if this is the first time through, just copy the results into the results.perm.list
					results.perm.list[[p]] <- one.perm
					}else{
					for(i in 1:length(one.perm)){
						results.perm.list[[p]][[i]] <- rbind(results.perm.list[[p]][[i]], one.perm[[i]])
						}
					}
				final.perm <- dim(results.perm.list[[p]][[1]])[1]
				j = j + 1
				} #end when we have enough permutations
			} #end looping over phenotypes


	# #plot the null distributions
	# layout.mat <- get.layout.mat(length(scanone.result), "landscape")
	# pdf(paste("Pair.Null.Dist.Top.", n.top.markers, ".Markers.", total.perm, ".Permutations.pdf", sep = ""), width = dim(layout.mat)[2]*4, height = dim(layout.mat)[1]*4)
	# layout(layout.mat)
	# for(i in 1:length(scanone.result)){
		# plot(density(results.perm.list[[i]][[1]][,"marker1:marker2"]/results.perm.list[[i]][[2]][,"marker1:marker2"]), main = paste("Null Dist. of Std Effects", names(scanone.result)[i], sep = "\n"), xlim = c(-5, 5))
		# }
	# dev.off()

	names(results.perm.list) <- names(scanone.result)
	data.obj$pairscan.perm <- results.perm.list
			
	return(data.obj) #and return it
	
}
