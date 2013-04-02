get.pairs.for.pairscan <-
function(geno, min.per.genotype, verbose = FALSE){
	
	if(verbose){
		cat("\nChecking marker pairs for genotype representation...\n")
		}
	
	
	check.linkage <- function(m1,m2,min.per.genotype) {
		t <- cbind.data.frame(as.factor(m1),as.factor(m2))
		colnames(t) <- c("m1","m2")
		reps <- table(t$m1,t$m2)
		too.few <- which(reps < min.per.genotype)
		if(length(too.few) >= 1) {
			return(FALSE) #pair failed check
			}else{
			return(TRUE) #pair passed check
			}
		}



	
	all.pairs <- pair.matrix(1:dim(geno)[2])
	all.pair.names <- pair.matrix(colnames(geno))
	
	check.one.pair <- function(pair){
		pass.checks <- check.linkage(m1 = geno[,pair[1]], m2 = geno[,pair[2]], min.per.genotype = min.per.genotype)
		return(pass.checks)
		}
	
	good.pairs <- apply(all.pairs, 1, check.one.pair)
	
	pairs.mat <- all.pair.names[which(good.pairs),]
	colnames(pairs.mat) <- c("marker1", "marker2")
	rownames(pairs.mat) <- NULL
	
	return(pairs.mat)
	
}
