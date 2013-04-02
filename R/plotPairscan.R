plotPairscan <-
function(data.obj, phenotype = "ET1", standardized = FALSE, pdf.label = "Pairscan.Regression.pdf"){
	
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
		
	num.pheno <- length(phenotype)
	
	#collect the results, so we can put them on the same scale
	all.results.mats <- list()
	min.x <- 0
	max.x <- 0
	#for each phenotype scanned
	for(p in 1:num.pheno){
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
				results.mat[marker1, marker2] <- as.numeric(pairscan.result[[p]][[1]][i,"marker1:marker2"])/as.numeric(pairscan.result[[p]][[2]][i,"marker1:marker2"])
				results.mat[marker2, marker1] <- as.numeric(pairscan.result[[p]][[1]][i,"marker1:marker2"])/as.numeric(pairscan.result[[p]][[2]][i,"marker1:marker2"])
				
			}else{
				results.mat[marker1, marker2] <- as.numeric(pairscan.result[[p]][[1]][i,"marker1:marker2"])
				results.mat[marker2, marker1] <- as.numeric(pairscan.result[[p]][[1]][i,"marker1:marker2"])
				}
			}
			
		results.mat[lower.tri(results.mat, diag = TRUE)] <- 0
		all.results.mats[[p]] <- results.mat
		min.x <- max(abs(c(max.x, results.mat)), na.rm = TRUE)*-1
		max.x <- max(abs(c(max.x, results.mat)), na.rm = TRUE)
		
		

		}	

		pdf(pdf.label)
		for(p in 1:num.pheno){
			myImagePlot(all.results.mats[[p]], xlab = "marker1", ylab = "marker2", main = phenotype[p], min.x = min.x, max.x = max.x)
			}
		dev.off()
		
		
	
	
}
