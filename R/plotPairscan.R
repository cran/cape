plotPairscan <-
function(data.obj, phenotype = NULL, standardized = TRUE, show.marker.labels = FALSE, show.chr = TRUE, label.chr = TRUE, pdf.label = "Pairscan.Regression.pdf", verbose = FALSE){
	
	#get the markers used in the pair scan and sort them.
	markers <- colnames(data.obj$geno.for.pairscan)
	marker.locale <- match(markers, colnames(data.obj$geno))
	sorted.markers <- markers[order(marker.locale)]

	markers.which <- match(sorted.markers, colnames(data.obj$geno))
	#get coordinates of the chromosome boundaries
	if(show.chr){
		chromosomes <- data.obj$chromosome[markers.which]
		u_chr <- unique(chromosomes)
		chr.boundaries <- apply(matrix(u_chr, ncol = 1), 1, function(x) max(which(chromosomes == x))) + 0.5
		chr.boundaries <- c(0, chr.boundaries)
		if(label.chr){
			chr.names <- unique(chromosomes)
			chr.names[which(chr.names == 0)] <- "c"
			}else{
			chr.names <- NULL	
			}
		}else{
		chr.boundaries <- NULL
		chr.names <- NULL
		}

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
	
	marker.pairs <- data.obj$pairscan.results[[1]][[1]][,1:2]

	
	marker1.pos <- match(marker.pairs[,1], sorted.markers)	
	marker2.pos <- match(marker.pairs[,2], sorted.markers)
	
	#collect the results, so we can put them on the same scale
	all.results.mats <- list()
	min.x <- 0
	max.x <- 0
	#for each phenotype scanned
	for(p in pheno.num){
		if(verbose){cat("\n", phenotype[p], ":\n", sep = "")}
		#build a results matrix
		results.mat <- matrix(0, length(markers), length(markers))
		colnames(results.mat) <- rownames(results.mat) <- sorted.markers
		
		if(standardized){
			pair.int <- as.numeric(pairscan.result[[p]][[1]][, "marker1:marker2"])/as.numeric(pairscan.result[[p]][[2]][,"marker1:marker2"])
			}else{
			pair.int <- as.numeric(pairscan.result[[p]][[1]][,"marker1:marker2"])
			}
		
		#and fill it in from the results in the table
		for(i in 1:length(pairscan.result[[p]][[1]][,1])){
			if(verbose){report.progress(i, length(pairscan.result[[p]][[1]][,1]))}
			#so we don't have to order the markers, put the effect
			#in the upper right and lower left. Then blank out the
			#upper right. Otherwise we get some entries in the top
			#triangle, and some in the bottom.
			results.mat[marker1.pos[i], marker2.pos[i]] <- pair.int[i]
			results.mat[marker2.pos[i], marker1.pos[i]] <- pair.int[i]	
			}
			
		# results.mat[lower.tri(results.mat, diag = TRUE)] <- 0
		# results.mat[upper.tri(results.mat, diag = TRUE)] <- NA
		diag(results.mat) <- NA
		all.results.mats[[p]] <- results.mat
		min.x <- max(abs(c(max.x, results.mat)), na.rm = TRUE)*-1
		max.x <- max(abs(c(max.x, results.mat)), na.rm = TRUE)
		
		}	
	


		pdf(pdf.label, width = 7, height = 6)
		for(p in 1:length(pheno.num)){
			myImagePlot(x = all.results.mats[[pheno.num[p]]], xlab = "marker1", ylab = "marker2", main = phenotype[p], min.x = min.x, max.x = max.x, show.labels = show.marker.labels, chromosome.coordinates = chr.boundaries, chr.names = chr.names)
			}
		dev.off()
		
		
	
	
}
