plotPairscan <-
function(data.obj, phenotype = NULL, standardized = FALSE, pdf.label = "Pairscan.Regression.pdf"){
	
	#get the markers used in the pair scan and sort them.
	markers <- colnames(data.obj$geno.for.pairscan)
	marker.locale <- match(markers, colnames(data.obj$geno))
	sorted.markers <- markers[order(marker.locale)]

	pairscan.result <- data.obj$pairscan.results
	
	if(is.null(pairscan.result)){
		stop("pairscan() must be run before plotPairscan()")
		}
	
	if(is.null(phenotype)){
		phenotype <- names(pairscan.result)
		}
		
	pheno.num <- which(names(pairscan.result) %in% phenotype)
	
	if(length(pheno.num) < length(phenotype)){
		not.found <- which(!(phenotype %in% names(pairscan.result)))
		message("I couldn't find the following phenotypes:")
		cat(phenotype[not.found], sep = "\n")
		stop()
		}
		
	num.pheno <- length(pheno.num)
	
	#collect the results, so we can put them on the same scale
	all.results.mats <- list()
	min.x <- 0
	max.x <- 0
	#for each phenotype scanned
	for(p in pheno.num){
		#build a results matrix
		results.mat <- matrix(0, length(markers), length(markers))
		colnames(results.mat) <- rownames(results.mat) <- sorted.markers
		#and fill it in from the results in the table
		for(i in 1:length(pairscan.result[[p]][[1]][,1])){
			marker1 <- pairscan.result[[p]][[1]][i,1]
			marker2 <- pairscan.result[[p]][[1]][i,2]
			
			#so we don't have to order the markers, put the effect
			#in the upper left and lower right. Then blank out the
			#lower right. Otherwise we get some entries in the top
			#triangle, and some in the bottom.
			if(standardized){
				# results.mat[as.character(marker1), as.character(marker2)] <- as.numeric(pairscan.result[[p]][[1]][i,"marker1:marker2"])/as.numeric(pairscan.result[[p]][[2]][i,"marker1:marker2"])
				# results.mat[as.character(marker2), as.character(marker1)] <- as.numeric(pairscan.result[[p]][[1]][i,"marker1:marker2"])/as.numeric(pairscan.result[[p]][[2]][i,"marker1:marker2"])
                results.mat[which(colnames(results.mat) == marker1), which(colnames(results.mat) == marker2)] <- as.numeric(pairscan.result[[p]][[1]][i, 
                  "marker1:marker2"])/as.numeric(pairscan.result[[p]][[2]][i, 
                  "marker1:marker2"])
                results.mat[which(colnames(results.mat) == marker2), which(colnames(results.mat) == marker1)] <- as.numeric(pairscan.result[[p]][[1]][i, 
                  "marker1:marker2"])/as.numeric(pairscan.result[[p]][[2]][i, 
                  "marker1:marker2"])
				
			}else{
				# results.mat[as.character(marker1), as.character(marker2)] <- as.numeric(pairscan.result[[p]][[1]][i,"marker1:marker2"])
				# results.mat[as.character(marker2), as.character(marker1)] <- as.numeric(pairscan.result[[p]][[1]][i,"marker1:marker2"])
                results.mat[which(colnames(results.mat) == marker1), which(colnames(results.mat) == marker2)] <- as.numeric(pairscan.result[[p]][[1]][i, 
                  "marker1:marker2"])
                results.mat[which(colnames(results.mat) == marker2), which(colnames(results.mat) == marker1)] <- as.numeric(pairscan.result[[p]][[1]][i, 
                  "marker1:marker2"])
				}
			}
			
		results.mat[lower.tri(results.mat, diag = TRUE)] <- 0
		marker.locale <- which(sorted.markers %in% colnames(data.obj$geno))
		colnames(results.mat) <- rownames(results.mat) <- data.obj$marker.names[marker.locale]
		all.results.mats[[p]] <- results.mat
		min.x <- max(abs(c(max.x, results.mat)), na.rm = TRUE)*-1
		max.x <- max(abs(c(max.x, results.mat)), na.rm = TRUE)
		
		}	

		pdf(pdf.label)
		for(p in 1:length(pheno.num)){
			myImagePlot(all.results.mats[[pheno.num[p]]], xlab = "marker1", ylab = "marker2", main = phenotype[p], min.x = min.x, max.x = max.x)
			}
		dev.off()
		
		
	
	
}
