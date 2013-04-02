genome.wide.threshold.1D <-
function(data.obj, n.perm = 1000, alpha.for.pairs = 0.05, alpha.for.covar = 0.01, scan.what = c("eigentraits", "raw.traits"), verbose = FALSE){
	
	require("evd")

	#calculate the numbers of markers, phenotypes and samples
	n.gene <- dim(data.obj$geno)[2]

	#pull out genotype and phenotype data for
	#clarity of code later.
	#If the user has not specified a scan.what,
	#from the outer function (singlescan.R),
	#default to eigentraits, basically, if eigen,
	#et, or ET are anywhere in the string, use the
	#eigentrais, otherwise, use raw phenotypes
	type.choice <- c(grep("eigen", scan.what), grep("ET", scan.what), grep("et", scan.what))
	if(length(type.choice) > 0){
		pheno <- data.obj$ET
		}else{
			pheno <- data.obj$pheno
			}

	gene <- data.obj$geno
	num.samples <- dim(pheno)[1]
	
	#create a matrix to hold permutation results
  	perm.max <- matrix(NA, ncol = length(pheno[1,]), nrow = n.perm)


      for (j in 1:n.perm) {
      	
      	if(verbose){
 	     	report.progress(current = j, total = n.perm, percent.text = 10, percent.dot = 2)
			}
		
		#shuffle the vector of individuals
        sampled.vector <- sample(1:num.samples)
		gene_perm <- gene[sampled.vector,]		
		
		#loop over traits and do all regressions on
		#markers at once
		
		get.stat <- function(regression){
			if(dim(summary(regression)$coefficients)[1] == 2){
				stat <- summary(regression)$coefficients[2,1]/summary(regression)$coefficients[2,2]
				}else{
					stat <- NA
					}
				return(stat)
			}
			
		stat.mat <- matrix(NA, nrow = n.gene, ncol = length(pheno[1,]))
  		for(et in 1:length(pheno[1,])){
  			regress.list <- apply(gene_perm, 2, function(x) lm(pheno[,et]~x))
            stat.mat[,et] <- as.vector(sapply(regress.list, get.stat))
  			}
                  
      		max.stat <- apply(stat.mat, 2, function(x) max(x, na.rm = TRUE))
            perm.max[j,] <- max.stat
    
        } #end permutations
        
	#apply the extreme value distribution to the results
	evd <- apply(perm.max, 2, function(x) fgev(x, std.err = FALSE))

	get.s <- function(evd.result, alpha){
		s <- qgev(1-alpha,loc=evd.result$estimate[1], scale=evd.result$estimate[2], shape=evd.result$estimate[3], lower.tail = TRUE)
		return(s)
		}
	
	
	s <- list()
	s[[1]] <- as.vector(sapply(evd, function(x) get.s(x, alpha.for.pairs)))
	s[[2]] <- as.vector(sapply(evd, function(x) get.s(x, alpha.for.covar)))
	
	#calculate one threshold over all phenotypes
	thresholds <- lapply(s, mean)

	data.obj$pairscan.thresh <- thresholds[[1]]
	data.obj$covar.thresh <- thresholds[[2]]
		
	if(verbose){
		cat("\n") #make sure the prompt is on the next line at the end of everything
		}
		
	return(data.obj)


}
