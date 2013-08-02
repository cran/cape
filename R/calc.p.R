calc.p <-
function(data.obj, pval.correction = c("holm", "fdr", "lfdr", "none")) {
	
	choice <- 0
	if(!is.null(data.obj$var.to.var.p.val)){
		message("\nIt appears that p-values have already been calculated for this data object.\n", sep = "")

		while(choice != "a" && choice != "b" && choice != "c"){
			choice <- readline(prompt = "Would you like to:\n\t(a) recalculate all p-values?\n\t(b) change the multiple testing correction?\n\t(c) cancel?\n")
			}
		if(choice == "c"){
			return(data.obj)
			}
		}
	
	# require("fdrtool")
	plot.null.transform = TRUE
	
	if(length(grep("h", pval.correction) > 0)){
		pval.correction <- "holm"
		}
		
	if(pval.correction != "holm" && pval.correction != "fdr" && pval.correction != "lfdr" && pval.correction != "none"){
		stop("pval.correction must be one of the following: 'holm', 'fdr', 'lfdr', 'none'")
		}
		
	
	if(choice == "a" || choice == 0){
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
	mat12.perm <- as.numeric(influences.perm[,3]) / as.numeric(influences.perm[,4])
	mat21.perm <- as.numeric(influences.perm[,5]) / as.numeric(influences.perm[,6])
	mat12.mat21.perm <- c(mat12.perm, mat21.perm)

	mat12 <- as.numeric(influences.org[,3]) / as.numeric(influences.org[,4])
	mat21 <- as.numeric(influences.org[,5]) / as.numeric(influences.org[,6])


	#I'm taking this correction out for now to test whether symmetrized genotypes get rid of the skewed distributions
	
	# #The m12/m21 distribution can often be really skewed, so here we implement
	# #a symmetrizing algorithm so we can calculate empirical p values
	# mat12.adj <- symmetrize.null(null.dist <- mat12.mat21.perm, emp.dist = mat12, plot.dist = plot.null.transform, plot.title = "Transformation.m12.Null.pdf")
	# mat21.adj <- symmetrize.null(null.dist <- mat12.mat21.perm, emp.dist = mat21, plot.dist = plot.null.transform, plot.title = "Transformation.m21.Null.pdf")	

	# adj.null <- mat12.adj$normalized.null
	# adj.m12 <- mat12.adj$normalized.empirical
	# adj.m21 <- mat21.adj$normalized.empirical


	adj.null <- mat12.mat21.perm
	adj.m12 <- mat12
	adj.m21 <- mat21


	#changed calculation of p value to account for the asymmetric m12/m21 distribution
	#I now calculate the p value based on above and below the median m12/m21
	#separately
	get.emp.p <- function(num.pair){
		emp.vals <- c(adj.m12[num.pair], adj.m21[num.pair])
		emp.p <- rep(NA, 2)
		for(e in 1:length(emp.vals)){
			if(emp.vals[e] < median(adj.null)){
				emp.p[e] <- length(which(adj.null <= emp.vals[e]))/length(which(adj.null <= median(adj.null)))
				}else{
				emp.p[e] <- length(which(adj.null >= emp.vals[e]))/length(which(adj.null >= median(adj.null)))
				}
			}
		return(emp.p)
		}


	
	# all.emp.p <- t(apply(matrix(1:n.pairs, ncol = 1), 1, function(x) get.emp.p(x)))
	all.emp.p <- t(apply(matrix(1:n.pairs, ncol = 1), 1, function(x) get.emp.p(x)))

	m12 <- matrix(c(marker.mat[,2],marker.mat[,1],as.numeric(as.matrix(influences.org[,3])),as.numeric(as.matrix(influences.org[,4])),(abs(as.numeric(influences.org[,3])) / as.numeric(influences.org[,4])),all.emp.p[,1]), ncol = 6)	
	colnames(m12) <- c("Source","Target","Effect","SE","|Effect|/SE","P_empirical")	

	m21 <- matrix(c(marker.mat[,1],marker.mat[,2],as.numeric(as.matrix(influences.org[,5])),as.numeric(as.matrix(influences.org[,6])),(abs(as.numeric(influences.org[,5])) / as.numeric(influences.org[,6])),all.emp.p[,2]), ncol = 6)
	colnames(m21) <- c("Source","Target","Effect","SE","|Effect|/SE","P_empirical")
	} #end if choice == "a"
	
	
	#adjust the p values
	if(choice == "b"){
		final.table <- data.obj$var.to.var.p.val
		final.table <- final.table[,-7]
		}else{
		final.table <- rbind(m12, m21)
		}
		
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
