direct.influence <-
function(data.obj, transform.to.phenospace = TRUE, verbose = FALSE) {


	data.obj$transform.to.phenospace <- transform.to.phenospace


	#calculate the direct influences for either the actual
	#tests or the permutations. Return a list with one element
	#for each phenotype with columns:
	#marker1, marker2, marker1.influence.coef, marker2.influence.coef, marker1.se, marker2.se
	
	direct.influence <- function(data.obj, perm){
		geno <- data.obj$geno.for.pairscan 
		n.gene <- dim(geno)[2]
		if(perm){
			scan.two.results <- data.obj$pairscan.perm		
			}else{
			scan.two.results <- data.obj$pairscan.results
			}
			
		n.pairs <- length(scan.two.results[[1]][[1]][,1]) 
		marker.mat <- scan.two.results[[1]][[1]][,1:2] 
		colnames(marker.mat) <- c("marker1", "marker2")
		orig.pheno <- data.obj$pheno
		data.obj$transform.to.phenospace <- transform.to.phenospace
	
		
		#if transform.to.phenospace is FALSE, figure out
		#if we are calculating the direct influences on
		#phenotypes or eigentraits
	
		if(!transform.to.phenospace){
			pheno.names <- names(data.obj$pairscan.results)	
			pheno.check <- match(pheno.names, colnames(data.obj$pheno))
			if(length(which(!is.na(pheno.check))) == 0){ #if we scanned eigentraits
				ET <- data.obj$ET					
				}else{
				ET <- data.obj$pheno	
				}					
			}else{
			ET <- data.obj$ET
			right.sing.vals <- data.obj$right.singular.vectors 
			diag.mat <- round(t(ET)%*%orig.pheno%*%right.sing.vals, 2) #get the singular value matrix
			pheno.names <- colnames(data.obj$pheno)			
			}			

			
		num.ET <- dim(ET)[2] 
		num.pheno <- length(pheno.names)
		coeff.names <- c("marker1", "marker2")
		
		#=========================================================
		# preallocate matrices to hold the statistics for each pair
		# The columns are for marker1, marker2, marker1.beta, 
		# marker2.beta, marker1.se and marker2.se (6 columns)
		# there is one of these for each phenotype
		#=========================================================
		stats.mat <- matrix(NA, ncol = 6, nrow = n.pairs)
		colnames(stats.mat) <- c("marker1", "marker2", "marker1.coef", "marker2.coef", "marker1.se", "marker2.se")
		stats.list <- vector(mode = "list", length = num.pheno)
		names(stats.list) <- pheno.names
		for(n in 1:num.pheno){
			stats.list[[n]] <- stats.mat
			}


		#This function grabs either the beta matrix or the se matrix
		#The result element determines which type of matrix will be 
		#returned: beta, result.element = 1, se, result.element = 2
		get.beta.mat <- function(scan.two.results, marker.pair.number, result.element){
			#the beta matrix is composed of the coefficients from each pairwise
			#marker model (except for the interaction coefficient)
			#we use all markers in each row and set the non-covariate entries to 0
			num.pheno <- length(scan.two.results)
			beta.mat <- matrix(NA, nrow = (dim(scan.two.results[[1]][[1]])[2]-3), ncol = num.pheno)
			for(ph in 1:num.pheno){
				beta.mat[,ph] <- as.numeric(scan.two.results[[ph]][[result.element]][marker.pair.number,3:(dim(scan.two.results[[ph]][[result.element]])[2]-1)])
				}				
			colnames(beta.mat) <-colnames(ET)
			rownames(beta.mat) <- colnames(scan.two.results[[ph]][[result.element]])[3:(dim(scan.two.results[[ph]][[result.element]])[2]-1)]
			return(beta.mat)	
			}
		



		#This function calculates a beta prime matrix from each
		#beta matrix by multiplying by the diagonal singular value 
		#matrix, and the transpose of the right singular values.	
		calculate.beta.prime <- function(beta.mat){
			#calculate the full beta prime matrix
			#if we are converting to phenotype space
			if(transform.to.phenospace){
				beta.prime <- beta.mat%*%diag.mat%*%t(right.sing.vals)
				#get only the stats for marker1 and marker2
				trunc.beta.prime <- matrix(beta.prime[(length(beta.prime[,1])-1):length(beta.prime[,1]),], nrow = 2)
				colnames(trunc.beta.prime) <- colnames(orig.pheno)
				}else{
				beta.prime <- beta.mat #otherwise, just use the straight beta matrix
				trunc.beta.prime <- matrix(beta.prime[(length(beta.prime[,1])-1):length(beta.prime[,1]),], nrow = 2)
				colnames(trunc.beta.prime) <- colnames(ET)
				}

			just.marker.beta <- vector(mode = "list", length = dim(trunc.beta.prime)[2])
			for(i in 1:length(just.marker.beta)){
				just.marker.beta[[i]] <- trunc.beta.prime[,i]
				}
			names(just.marker.beta) <- colnames(trunc.beta.prime)
			return(just.marker.beta)
			}



		#This function calculates a beta se prime matrix from each
		#beta matrix by multiplying by the beta se, squared diagonal singular value 
		#matrix, and the square of the transpose of the right singular values.	
		calculate.se.prime <- function(se.mat){

			sqrd.se.mat <- se.mat ^ 2
			if(transform.to.phenospace){
				temp.mat <- diag.mat %*% t(right.sing.vals)
				sqrd.temp.mat <- (temp.mat ^  2)
				se.prime <- sqrd.se.mat %*% sqrd.temp.mat
				colnames(se.prime) <- colnames(orig.pheno)
				}else{
					se.prime <- sqrd.se.mat
					colnames(se.prime) <- colnames(ET)
					}

			se.prime <- sqrt(se.prime)

			just.marker.se <- vector(mode = "list", length = dim(se.prime)[2])
			for(i in 1:length(just.marker.se)){
				just.marker.se[[i]] <- se.prime[((dim(se.prime)[1]-1):dim(se.prime)[1]),i]
				}
			names(just.marker.se) <- colnames(se.prime)
		
			return(just.marker.se)
		}



		#=========================================================
		# Calculation of direct influences
		# 1. Calculation of beta prime
		#=========================================================

		#get the beta matrix for each pair of markers
		beta.matrices <- vector(mode = "list", length = n.pairs)
		for(p in 1:length(marker.mat[,1])){
			beta.matrices[[p]] <- get.beta.mat(scan.two.results, p, 1)
			}
	
		#calculate beta prime for each of the beta matrices
		all.beta.prime <- lapply(beta.matrices, calculate.beta.prime)
		
		#add these into the final stats.matrix
		for(i in 1:length(stats.list)){
			stats.list[[i]][,1:2] <- marker.mat
			stats.list[[i]][,3:4] <- t(sapply(all.beta.prime, function(x) x[[i]]))
			}

	
		#=========================================================
		# Calculation of direct influences
		# 2. Calculation of beta prime's std error using error prop
		#=========================================================

		#get the beta matrix for each pair of markers
		se.matrices <- vector(mode = "list", length = n.pairs)
		for(p in 1:length(marker.mat[,1])){
			se.matrices[[p]] <- get.beta.mat(scan.two.results, p, 2)
			}
	
		#calculate beta prime for each of the beta matrices
		all.se.prime <- lapply(se.matrices, calculate.se.prime)
		
		#add these into the final stats.matrix
		for(i in 1:length(stats.list)){
			stats.list[[i]][,5:6] <- t(sapply(all.se.prime, function(x) x[[i]]))
			}

		#=========================================================
		# Generate final results list object
		#=========================================================

		if(perm){
			data.obj$var.to.pheno.influence.perm <- stats.list			
			}else{
			data.obj$var.to.pheno.influence <- stats.list				
				}

			return(data.obj)
			}



	

	#separate out the markers from each pair
	#for each marker in each context, give it an 
	#influence coefficient, influence se, t stat
	#and |t stat|
	
	direct.influence.stat <- function(data.obj, perm){
		
		if(perm){
			marker.stats <- data.obj$var.to.pheno.influence.perm	
			}else{
			marker.stats <- data.obj$var.to.pheno.influence
			}

		markers <- colnames(data.obj$geno.for.pairscan)
		pheno.names <- names(marker.stats)
		stats.list <- vector(length = length(pheno.names), mode = "list")
		names(stats.list) <- pheno.names
		
		marker.incidence <- apply(matrix(markers, ncol = 1), 1, function(x) length(which(marker.stats[[1]][,1:2] == x)))
		#preallocate a matrix that will hold the statistics for each marker
		#in each pair context. The columns are marker, coef, se, t.stat
		#|t.stat|, emp.p (6 columns)
		if(perm){
			stats.mat <- matrix(NA, nrow = sum(marker.incidence), ncol = 5)
			colnames(stats.mat) <- c("marker", "coef", "se", "t.stat", "|t.stat|")
			}else{
			stats.mat <- matrix(NA, nrow = sum(marker.incidence), ncol = 6)
			colnames(stats.mat) <- c("marker", "coef", "se", "t.stat", "|t.stat|", "emp.p")			
			}
		
		stats.list <- vector(mode = "list", length = length(pheno.names))
		names(stats.list) <- pheno.names
		for(i in 1:length(stats.list)){
			stats.list[[i]] <- stats.mat
			}
		
		#for each phenotype, go through the marker names and
		#collect the statistics from var.to.pheno.influence 
		#for when the marker was in position 1 and in position 2
		for(ph in 1:length(pheno.names)){
			
			#go through each marker and take out its statistics
			for(m in markers){

				#find out where the marker was in position 1
				m1.locale <- which(marker.stats[[ph]][,1] == m)
				
				if(length(m1.locale) > 0){
					beta.coef1 <- as.numeric(marker.stats[[ph]][m1.locale, "marker1.coef"])
					se1 <- as.numeric(marker.stats[[ph]][m1.locale, "marker1.se"])
					t.stat1 <- beta.coef1/se1
					stat.section1 <- matrix(c(rep(m, length(beta.coef1)), beta.coef1, se1, t.stat1, abs(t.stat1)), ncol = 5)

					start.row <- min(which(is.na(stats.list[[ph]][,1])))
					stats.list[[ph]][start.row:(start.row + dim(stat.section1)[1] - 1),1:dim(stat.section1)[2]] <- stat.section1
					}

				#find out where the marker was in position 2
				m2.locale <- which(marker.stats[[ph]][,2] == m)

				if(length(m2.locale) > 0){
					beta.coef2 <- as.numeric(marker.stats[[ph]][m2.locale, "marker2.coef"])
					se2 <- as.numeric(marker.stats[[ph]][m2.locale, "marker2.se"])
					t.stat2 <- beta.coef2/se2
					stat.section2 <- matrix(c(rep(m, length(beta.coef2)), beta.coef2, se2, t.stat2, abs(t.stat2)), ncol = 5)
				
					start.row <- min(which(is.na(stats.list[[ph]][,1])))
					stats.list[[ph]][start.row:(start.row + dim(stat.section2)[1] - 1),1:dim(stat.section2)[2]] <- stat.section2
					}
	
				}

			}
		
		if(perm){
			data.obj$var.to.pheno.test.stat.perm <- stats.list
			}else{
			data.obj$var.to.pheno.test.stat <- stats.list
			}
			
		return(data.obj)
	}
	
	
	
	direct.influence.epcal<- function(data.obj){
	
		stat <- data.obj$var.to.pheno.test.stat
		perm.test.stat <- data.obj$var.to.pheno.test.stat.perm
		
		get.p <- function(val, dist){
			pval <- length(which(dist >= val))/length(dist)
			return(pval)
			}

		#caclualte an empirical p value for each entry in stat
		for(ph in 1:length(stat)){
			comp.dist <- as.numeric(perm.test.stat[[ph]][,"|t.stat|"]) #compare each calculated value to the permuted distribution
			stat[[ph]][,"emp.p"] <- apply(matrix(as.numeric(stat[[ph]][,"|t.stat|"]), ncol = 1), 1, function(x) get.p(x, comp.dist))
			}
		
		#replace the direct influence tables with the tables with the added empirical p values
		data.obj$var.to.pheno.test.stat <- stat
		return(data.obj)
		}


	
	direct.influence.ep.adj <- function(data.obj){

		stat <- data.obj$var.to.pheno.test.stat
		
		for(ph in 1:length(stat)){
			emp.pvals <- as.numeric(stat[[ph]][,"emp.p"])
			adj.p <- p.adjust(emp.pvals)
			new.table <- cbind(stat[[ph]], adj.p)
			ordered.table <- new.table[order(as.numeric(new.table[,"|t.stat|"]), decreasing = TRUE),] #order by increasing t statistic
			data.obj$var.to.pheno.test.stat[[ph]] <- ordered.table
			}

		return(data.obj)
		}



	direct.influence.max <- function(data.obj){

		stat <- data.obj$var.to.pheno.test.stat

		markers <- unique(stat[[1]][,1])
		max.stat.list <- vector(length(stat), mode = "list")
		names(max.stat.list) <- names(stat)
		
		#go through each of the phenotypes and find the maximum influence of each marker
		for(ph in 1:length(stat)){
			max.stat.table <- NULL
			for(m in markers){
				marker.locale <- which(stat[[ph]][,1] == m)
				temp.table <- matrix(stat[[ph]][marker.locale,], nrow = length(marker.locale))
				colnames(temp.table) <- colnames(stat[[ph]])
				max.stat.locale <- which(as.numeric(temp.table[,"|t.stat|"]) == max(as.numeric(temp.table[,"|t.stat|"])))
				max.stats <- matrix(temp.table[max.stat.locale,], nrow = length(max.stat.locale))
				max.stat.table <- rbind(max.stat.table, max.stats[1,]) #if there is more than one max, just take one
				colnames(max.stat.table) <- colnames(stat[[1]])
				}		
			max.stat.list[[ph]] <- max.stat.table
			}
	
		data.obj$max.var.to.pheno.influence <- max.stat.list
		return(data.obj)
		}

	
	direct.influence.ep.adj.post.max <- function(data.obj){

		stat <- data.obj$max.var.to.pheno.influence
		
		for(ph in 1:length(stat)){
			emp.pvals <- as.numeric(stat[[ph]][,"emp.p"])
			adj.p <- p.adjust(emp.pvals)
			new.table <- cbind(stat[[ph]], adj.p)
			ordered.table <- new.table[order(as.numeric(new.table[,"|t.stat|"]), decreasing = TRUE),] #order by increasing t statistic
			data.obj$max.var.to.pheno.influence[[ph]] <- ordered.table
			}

		return(data.obj)
		}

	if(verbose){cat("calculating direct influence of variants...\n")}
	data.obj <- direct.influence(data.obj, perm = FALSE)
	if(verbose){cat("calculating direct influence of permutations...\n")}
	data.obj <- direct.influence(data.obj, perm = TRUE)
	if(verbose){cat("calculating p values from permutations...\n")}
	data.obj <- direct.influence.stat(data.obj, perm = FALSE)
	data.obj <- direct.influence.stat(data.obj, perm = TRUE)
	data.obj <- direct.influence.epcal(data.obj)
	if(verbose){cat("adjusting p values with Holm's stepdown procedure...\n")}	
	# data.obj <- direct.influence.ep.adj(data.obj)	
	data.obj <- direct.influence.max(data.obj)
	data.obj <- direct.influence.ep.adj.post.max(data.obj)
	
	
	return(data.obj)
}
