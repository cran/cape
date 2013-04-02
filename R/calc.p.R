calc.p <-
function(data.obj) {
	
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


	#### Use when combining across permutations#####
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

	m12 <- cbind(marker.mat[,2],marker.mat[,1],as.numeric(as.matrix(influences.org[,3])),as.numeric(as.matrix(influences.org[,4])),(abs(as.numeric(influences.org[,3])) / as.numeric(influences.org[,4])),all.emp.p[,1])	
	colnames(m12) <- c("Source","Target","Effect","SE","|Effect|/SE","P_empirical")
	
	#adjust the p values
	p.adjusted <- p.adjust(as.numeric(m12[,"P_empirical"]), method = "holm")
	m12 <- cbind(m12, p.adjusted)


	m21 <- cbind(marker.mat[,1],marker.mat[,2],as.numeric(as.matrix(influences.org[,5])),as.numeric(as.matrix(influences.org[,6])),(abs(as.numeric(influences.org[,5])) / as.numeric(influences.org[,6])),all.emp.p[,2])	
	colnames(m21) <- c("Source","Target","Effect","SE","|Effect|/SE","P_empirical")
	p.adjusted <- p.adjust(as.numeric(m21[,"P_empirical"]), method = "holm")
	m21 <- cbind(m21, p.adjusted)


	final.table <- rbind(m12, m21)
	final.table <- final.table[order(final.table[,"|Effect|/SE"], decreasing = TRUE),]

	data.obj$var.to.var.p.val <- final.table
	
	return(data.obj)
}
