get.linearly.independent <-
function(data.obj){


	matrixX <- data.obj$geno.for.pairscan

	if(dim(matrixX)[2] == 1){
		return(matrixX)
		}
	
	#use precision to the 3rd decimal place
	orig.matrixX <- matrixX
	matrixX <- round(matrixX, 3)

	#find the markers without variation
	num.geno <- apply(matrixX, 2, function(x) length(unique(x[!is.na(x)])))
	rejected.markers <- names(num.geno[num.geno < 2])
	good.markers <- num.geno[num.geno >= 2]

	all.pairs <- pair.matrix(names(good.markers))

	get.cor <- function(pair){
		geno1 <- matrixX[,pair[1]]; geno2 <- matrixX[,pair[2]]
		good.vals <- sort(unique(intersect(which(!is.na(geno1)), which(!is.na(geno2)))))
		if(length(good.vals) > 3){
			if(var(matrixX[good.vals,pair[1]]) > 0 && var(matrixX[good.vals,pair[2]]) > 0){
				pair.cor <- cor(matrixX[,pair[1]], matrixX[,pair[2]], use = "complete.obs")
				}else{
					pair.cor <- NA
					}
			}else{
			pair.cor <- NA
			}
		return(pair.cor)
		}

	all.cor <- apply(all.pairs, 1, get.cor)
	
	perfect.cor <- which(abs(all.cor) == 1)
	#if there are markers with perfect correlation, 
	#remove the first one of the pair
	if(length(perfect.cor) > 0){
		for(i in perfect.cor){
			rejected.markers <- c(rejected.markers, all.pairs[i,1])
			}
		}

	
	rejected.markers <- as.vector(rejected.markers)
	if(length(rejected.markers) > 0){
		rej.markers <- match( sort(unique(rejected.markers)), colnames(matrixX))
		final.geno <- matrixX[,-sort(unique(rej.markers))]
		}else{
			rej.markers <- NULL
			final.geno <- matrixX
			}

	results <- list(final.geno, rej.markers)
	names(results) <- c("independent.markers", "rejected.markers")
	return(results)	

}
