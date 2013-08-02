select.markers.for.pairscan <-
function(data.obj, use.pairs.threshold = FALSE, pairscan.thresh = NULL, specific.markers = NULL, num.markers = NULL, start.thresh = 4, tolerance = 10, verbose = TRUE){
	
	# if(!is.null(num.markers)){
		# use.pairs.threshold = TRUE
		# }

	if(!use.pairs.threshold){
		pairscan.thresh <- NULL
		}

	#take our various parts of the data object
	#for clarity of code later on
	scanone.result <- data.obj$singlescan.results
	covar.flags <- data.obj$covar.flags
	
	#we need to use the covariates determined from the 1D scan
	#in the 2D scan, so stop if the 1D scan has not been performed
	if(length(scanone.result) == 0){
		stop("singlescan() must be run before selecting markers for pairscan\n")
		}	

	orig.geno <- data.obj$geno
	orig.covar.flags <- data.obj$covar.flags
	

	locus.names <- colnames(orig.geno)
	phenotype.names <- names(scanone.result)
	

	#if we are filtering by the pairs threshold, we 
	#need to find all alleles that exceed the significance
	#threshold for any of the phenotypes. All these alleles
	#get put into a 2D matrix. The rest of the algorithm is
	#identical to cape
	#if we are not thresholding, just flatten the entire
	#array to 2D
			
		
		#============================================================================================
		# Some functions to help parse data structures
		#============================================================================================
		get.sig.allele <- function(result.mat, pairs.thresh){
			sig.alleles <- which(abs(result.mat[,"t.stat"]) >= pairs.thresh, arr.ind = TRUE)
			return(sig.alleles)
			}
		
		# get.allele.row.pos <- function(allele.row){
			# locus.locale <- which(locus.names %in% allele.row[1])
			# allele.locale <- which(allele.names %in% allele.row[2])
			# return(c(locus.locale, allele.locale))
			# }
		
		get.named.alleles <- function(specific.markers){
			marker.locale <- lapply(scanone.result, function(x) which(specific.markers %in% rownames(x)))
			for(i in 1:length(marker.locale)){
				names(marker.locale[[i]]) <- specific.markers
				}
			return(marker.locale)
			}
		
		get.alleles <- function(above.thresh, locus){
			alleles <- c(sort(unique(unlist(lapply(above.thresh, function(x) x[which(x[,"locus"] == locus), "allele"])))))
			allele.list <- as.vector(alleles, mode = "list")
			names(allele.list) <- rep(locus, length(alleles))
			return(allele.list)
			}

		# extract.allele.geno <- function(allele.row){
			# locus <- orig.geno[,,as.numeric(allele.row[1])]
			# just.allele <- matrix(locus[,as.numeric(allele.row[2])], ncol = 1)
			# allele.label <- paste(locus.names[as.numeric(allele.row[1])], allele.names[as.numeric(allele.row[2])], sep = "_")	
			
			# just.covar <- orig.covar.flags[as.numeric(allele.row[1]),,]
			# just.covar.allele <- just.covar[,as.numeric(allele.row[2])]
			
			# result <- list(allele.label, just.allele, just.covar.allele)
			# names(result) <- c("allele.label", "genotype", "covar")
			# return(result)
			# }
					
			
		get.unique.sig <- function(above.thresh){
			all.loci.names <- sort(unique(unlist(sapply(above.thresh, function(x) names(x)))))
			locus.locale <- which(colnames(orig.geno) %in% all.loci.names)
			results <- list(orig.geno[,locus.locale], orig.covar.flags[locus.locale,])
			names(results) <- c("geno.for.pairscan", "covar.for.pairscan")
			return(results)
			}
								
		#one big function for selecting markers. 
		#this is repeated multiple times if num.markers is specified
		
		get.markers <- function(pairscan.thresh){
			#if no markers have been specified, get above.thresh using the pairs.threshold
			if(is.null(specific.markers)){ 
				# cat("specific.markers is NULL\n")
				if(is.null(pairscan.thresh) || use.pairs.threshold == FALSE){
					pairscan.thresh <- 0
					}
					if(is.null(pairscan.thresh)){
						#start at the minimum maximum value for all results
						pairs.thresh <- min(sapply(scanone.result, function(x) max(x, na.rm = TRUE)))
						}else{
						pairs.thresh <- pairscan.thresh
						}
					data.obj$pairscan.thresh <- pairs.thresh
					above.thresh <- lapply(scanone.result, function(x) get.sig.allele(x, pairs.thresh))
					if(sum(sapply(above.thresh, length)) == 0){stop("There are no markers above this threshold. Please choose a lower starting threshold.")}
					#make sure there is at least one instance of each non-allelic covariate
					covar.names <- rownames(covar.flags)[unique(apply(covar.flags, 2, function(x) which(x == 1)))]
					if(length(covar.names) > 0){
						for(i in 1:length(covar.names)){
							covar.locale <- lapply(above.thresh, function(x) which(names(x) %in% covar.names[i]))
							missing.covar <- lapply(covar.locale, function(x) if(length(x) == 0){TRUE}else{FALSE})
							for(j in 1:length(covar.locale)){
								if(missing.covar[[j]]){ #if there is no covariate in the section, find it and add it
									covar.pos <- which(rownames(covar.flags) == covar.names[i])
									above.thresh[[j]] <- c(above.thresh[[j]], covar.pos)
									names(above.thresh[[j]])[length(above.thresh[[j]])] <- covar.names[i]
									}
								}
							}
						}

					
					}else{ 
					#if marker names have been supplied, 
					#create above.thresh using these
					above.thresh <- get.named.alleles(specific.markers)
					}
					
				#if there are non-allelic covariates specified, make sure there
				#is only one instance of each.
				if(length(covar.names) > 0){
					for(i in 1:length(covar.names)){
						covar.locale <- lapply(above.thresh, function(x) which(names(x) %in% covar.names[i]))
						for(j in 1:length(covar.locale)){
							if(length(covar.locale[[j]]) > 1){
								above.thresh[[j]] <- above.thresh[[j]][-covar.locale[[j]][2:length(covar.locale[[j]])]]
								}
							}
						}
					}
				
				matrices.for.pairscan <- get.unique.sig(above.thresh)
				pair.geno <- matrices.for.pairscan$geno.for.pairscan
				pair.covar.flags <- matrices.for.pairscan$covar.for.pairscan
	
	
				data.obj$geno.for.pairscan <- pair.geno
				data.obj$covar.for.pairscan <- pair.covar.flags	
			
			
				#make sure that the genotype matrix is linearly independent
				geno.ind <- get.linearly.independent(data.obj)
				return(list(geno.ind[[1]], geno.ind[[2]], data.obj))
				}
			
		#============================================================================================
		# end functions			
		#============================================================================================
			
			
			
			if(is.null(num.markers)){
				test.obj <- get.markers(pairscan.thresh)
				data.obj <- test.obj[[3]]
				geno <- test.obj[[1]]
				rejected.markers <- test.obj[[2]]
				num.markers.selected <- dim(geno)[2]
				pair.geno <- data.obj$geno.for.pairscan
				pair.covar.flags <- data.obj$covar.for.pairscan
				if(!is.null(pairscan.thresh)){
					data.obj$pairscan.thresh <- data.obj$pairscan.thresh
					data.obj$alpha.for.pairs <- paste("t =", pairscan.thresh)
					}
				}else{
				if(verbose){cat("Finding a threshold...\n")}
				if(verbose){cat("thresh\tnum.markers\n")}
				#start with a mid-level threshold
				pairscan.thresh <- start.thresh
				update.val <- 0.01
				less.than <- 1; greater.than <- 1
				iteration <- 1
				test.obj <- get.markers(pairscan.thresh)
				pairscan.thresh <- test.obj[[3]]$pairscan.thresh
				num.markers.selected <- dim(test.obj[[1]])[2]
				if(verbose){cat(pairscan.thresh, "\t", num.markers.selected, "\n")}
				while((num.markers.selected < num.markers - tolerance) || num.markers.selected > num.markers + tolerance){
					if(num.markers.selected < num.markers){
						pairscan.thresh <- pairscan.thresh - update.val
						less.than <- less.than + 1
						}
					if(num.markers.selected > num.markers){
						pairscan.thresh <- pairscan.thresh + update.val
						greater.than <- greater.than + 1
						}
					test.obj <- get.markers(pairscan.thresh)
					num.markers.selected <- dim(test.obj[[1]])[2]
					if(verbose){cat(pairscan.thresh, "  ", num.markers.selected, "\n")}
					iteration <- iteration + 1
					if(greater.than > 2 & less.than > 2){
						stop("I can't get a number of markers within the tolerance. \nPlease adjust the tolerance or the number of markers desired")
						}
					}
				data.obj <- test.obj[[3]]
				geno <- test.obj[[1]]
				rejected.markers <- test.obj[[2]]
				num.markers.selected <- dim(geno)[2]
				pair.geno <- data.obj$geno.for.pairscan
				pair.covar.flags <- data.obj$covar.for.pairscan
				data.obj$alpha.for.pairs <- paste("t =", pairscan.thresh)
				cat("\n")
				}
			


			if(length(rejected.markers) > 0){
				pair.covar.flags <- pair.covar.flags[match(colnames(geno), rownames(pair.covar.flags)),]
				message("\n", length(rejected.markers), " marker(s) rejected due to linear non-independence.\n For more information see markers.removed.txt")
				write.table(colnames(pair.geno)[sort(rejected.markers)], "markers.removed.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
				}
						
		
			#replace the geno and covar flags objects in data.obj with the 
			#linear independent versions
			#sort the genotype matrix and covariate flags before
			#returning them.
			
			marker.order <- match(colnames(geno), colnames(pair.geno))
			geno <- geno[,order(marker.order)]

			marker.order <- match(rownames(pair.covar.flags), colnames(pair.geno))
			covar.flags <- pair.covar.flags[order(marker.order),]
						
			
			data.obj$geno.for.pairscan <- geno
			data.obj$covar.for.pairscan <- covar.flags

			
			cat(dim(geno)[2], "markers were selected for the pairscan.\n")
			cat("This makes", choose(dim(geno)[2], 2), "possible pairs.\n")

		return(data.obj)
	
}
