calc.p <-
function(data.obj, pval.correction = c("holm", "fdr", "lfdr", "none")) {
	
	# require("fdrtool")
	
	if(length(grep("h", pval.correction) > 0)){
		pval.correction <- "holm"
		}
		
	if(pval.correction != "holm" && pval.correction != "fdr" && pval.correction != "lfdr" && pval.correction != "none"){
		stop("pval.correction must be one of the following: 'holm', 'fdr', 'lfdr', 'none'")
		}
		
	
	influences.org <- data.obj$var.to.var.influences
	influences.perm <- data.obj$var.to.var.influences.perm
		
	if(is.null(influences.org)){
		stop("error.prop() with perm = FALSE must be run before running calc.p()")
		}
		
	if(is.null(influences.perm)){
		stop("error.prop() with perm = TRUE must be run before running calc.p()")
		}
		



	n.gene <- dim(data.obj$geno.for.pairscan)[2] #get the number of genes used in the pair scan
	n.pairs <- dim(data.obj$pairscan.results[[1]][[1]])[1] #the number of pairs scanned in the pairscan
    
    	
    marker.mat <- influences.org[,1:2] #a matrix listing the names of all marker combinations
    colnames(marker.mat) <- c("marker1", "marker2")


	#### Combinine across permutations#####
	#get the t statistics for all permutations
	mat12.perm <- abs(as.numeric(influences.perm[,3])) / as.numeric(influences.perm[,4])
	mat21.perm <- abs(as.numeric(influences.perm[,5])) / as.numeric(influences.perm[,6])
	mat12.mat21.perm <- c(mat12.perm, mat21.perm)

	mat12 <- abs(as.numeric(influences.org[,3])) / as.numeric(influences.org[,4])
	mat21 <- abs(as.numeric(influences.org[,5])) / as.numeric(influences.org[,6])


	get.emp.p <- function(num.pair){
		m12.emp.p <- length(which(mat12.mat21.perm >= mat12[num.pair]))/length(mat12.mat21.perm)
		m21.emp.p <- length(which(mat12.mat21.perm >= mat21[num.pair]))/length(mat12.mat21.perm)
		return(c(m12.emp.p, m21.emp.p))
		}
	
	all.emp.p <- t(apply(matrix(1:n.pairs, ncol = 1), 1, function(x) get.emp.p(x)))

	m12 <- matrix(c(marker.mat[,2],marker.mat[,1],as.numeric(as.matrix(influences.org[,3])),as.numeric(as.matrix(influences.org[,4])),(abs(as.numeric(influences.org[,3])) / as.numeric(influences.org[,4])),all.emp.p[,1]), ncol = 6)	
	colnames(m12) <- c("Source","Target","Effect","SE","|Effect|/SE","P_empirical")	

	m21 <- matrix(c(marker.mat[,1],marker.mat[,2],as.numeric(as.matrix(influences.org[,5])),as.numeric(as.matrix(influences.org[,6])),(abs(as.numeric(influences.org[,5])) / as.numeric(influences.org[,6])),all.emp.p[,2]), ncol = 6)
	colnames(m21) <- c("Source","Target","Effect","SE","|Effect|/SE","P_empirical")
	

	#adjust the p values
	final.table <- rbind(m12, m21)
	if(pval.correction == "none"){
		p.adjusted <- as.numeric(final.table[,"P_empirical"])
		final.table <- cbind(final.table, p.adjusted)
		}
	if(pval.correction == "holm"){
		p.adjusted <- p.adjust(as.numeric(final.table[,"P_empirical"]), method = "holm")
		final.table <- cbind(final.table, p.adjusted)
		}
		
	if(pval.correction == "fdr" || pval.correction == "lfdr"){
		fdr.out <- fdrtool(as.numeric(final.table[,"P_empirical"]), statistic = "pvalue", plot = FALSE, verbose = FALSE, cutoff.method = "fndr")
		if(pval.correction == "lfdr"){
			lfdr <- fdr.out$lfdr
			final.table <- cbind(final.table, lfdr)
			}else{
			qval <- fdr.out$qval
			final.table <- cbind(final.table, qval)	
			}
		}


	final.table <- final.table[order(final.table[,"|Effect|/SE"], decreasing = TRUE),]

	data.obj$var.to.var.p.val <- final.table
	
	return(data.obj)
}
