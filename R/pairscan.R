pairscan <-
function(data.obj, scan.what = c("eigentraits", "raw.traits"), n.perm = NULL, min.per.genotype = NULL, max.pair.cor = NULL, verbose = FALSE, num.pairs.limit = 1e4, num.perm.limit = 1e7) {

	if(is.null(n.perm)){
		stop("The number of permutations must be specified.")
		}

	#take our various parts of the data object
	#for clarity of code later on.
	#taking the scanone.results just tells us which
	#phenotypes were scanned
	scanone.result <- data.obj$singlescan.results
	
	#we need to use the covariates determined from the 1D scan
	#in the 2D scan, so stop if the 1D scan has not been performed
	if(length(scanone.result) == 0){
		stop("singlescan() must be run on this object before pairscan()\n")
		}	


	#If the user does not specify a scan.what, 
	#default to eigentraits, basically, if eigen,
	#et, or ET are anywhere in the string, use the
	#eigentraits, otherwise, use raw phenotypes
	type.choice <- c(grep("eig", scan.what), grep("ET", scan.what), grep("et", scan.what)) #look for any version of eigen or eigentrait, the user might use.
	if(length(type.choice) > 0){ #if we find any, use the eigentrait matrix
		pheno <- data.obj$ET
		if(is.null(pheno)){stop("There are no eigentraits. Please set scan.what to raw.traits, or run get.eigentraits().")}
		}else{
			pheno <- data.obj$pheno #otherwise, use the raw phenotype matrix
			}

	geno <- data.obj$geno.for.pairscan
	covar.flags <- data.obj$covar.for.pairscan

	if(is.null(geno)){
		stop("select.markers.for.pairscan() must be run before pairscan()")
		}
		
		
	num.markers <- dim(geno)[2]
	
	#fill in a matrix to index the marker pairs
	marker.matrix <- pair.matrix(colnames(geno))
	pared.marker.mat <- get.pairs.for.pairscan(geno, min.per.genotype = min.per.genotype, max.pair.cor = max.pair.cor, verbose = verbose)

	num.pairs <- dim(pared.marker.mat)[1]
	
	if(num.pairs == 0){
		stop("There are no pairs to test. Try lowering min.per.genotype or raising max.pair.cor.")
		}

	if(!is.null(num.pairs.limit) && num.pairs > num.pairs.limit){
		cat("\nThe number of pairs (",num.pairs,") exceeds ", num.pairs.limit, ".\n", sep = "")
		go.on <- readline(prompt = "Do you want to continue (y/n)?\n")
		if(length(grep("n", go.on))){
			message("Stopping pairwise scan...\n")
			return(data.obj)
		}else{
			cat("Continuing pairwise scan...\n")
		}
	}

	
	
	#make a list to hold the results. 
	#Tables from each of the phenotypes will be
	#put into this list
	results.list <- list()
	results.perm.list <- list()


	#run one.pairscan for each phenotype with results in scanone.result
	# if n.perm > 1 then call one.pairscan.perm()
	for(p in 1:length(scanone.result)){ 

		if(verbose){
			cat("\nScanning phenotype ", colnames(pheno)[p], ":\n", sep = "")
			}
				
		pairscan.results <- one.pairscan(phenotype.vector = pheno[,p], genotype.matrix = geno, covar.vector = covar.flags[,p], pairs.matrix = pared.marker.mat, n.perm = 0, verbose = verbose)
		results.list[[p]] <- pairscan.results[[1]]
		# results.perm.list[[p]] <- pairscan.results[[2]]
		} #end looping over phenotypes
	 
	 
	#generate the null distribution using pairscan.null
	num.pairs <- choose(dim(data.obj$geno.for.pairscan)[2], 2)
	total.perm = n.perm*num.pairs

	#calculate the number of top markers to use based
	#on the thresholding of the singlescan
	n.top.markers <- dim(data.obj$geno.for.pairscan)[2]
	
	data.obj <- pairscan.null(data.obj, scan.what = scan.what, total.perm = total.perm, n.top.markers = n.top.markers, verbose = verbose)

	names(results.list) <- names(scanone.result)


	data.obj$"pairscan.results" <- results.list	  #add the results to the data object



	return(data.obj) #and return it
	
}
