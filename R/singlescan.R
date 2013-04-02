singlescan <-
function(data.obj, n.perm = NULL, covar = NULL, scan.what = c("eigentraits", "raw.traits"), auto.covar.selection = FALSE, alpha.for.covar = 0.01, alpha.for.pairs = 0.1, verbose = FALSE) {
	
	if(is.null(n.perm)){
		stop("The number of permutations must be specified.")
		}
	
	
	if(n.perm < 2){
		stop("The number of permutations must be at least 2.")
		}
		
		
	#If the user does not specify a scan.what, 
	#default to eigentraits, basically, if eigen,
	#et, or ET are anywhere in the string, use the
	#eigentrais, otherwise, use raw phenotypes
	type.choice <- c(grep("eigen", scan.what), grep("ET", scan.what), grep("et", scan.what))
	if(length(type.choice) > 0){
		pheno <- data.obj$ET
		if(is.null(pheno)){
			message("\nI am supposed to scan eigentraits, but if looks as if get.eigentraits() has not been run yet.")
			}
		}else{
			pheno <- data.obj$pheno
			}


	gene <- data.obj$geno
	
	#we do not scan markers on the sex chromosomes
	#take these out here.
	x.locale <- grep("X", data.obj$chromosome, ignore.case = TRUE)
	if(length(x.locale) > 0){
		message("\nRemoving markers on the X chromosome")
		gene <- gene[,-x.locale]
		data.obj$chromosome <- data.obj$chromosome[-x.locale]
		data.obj$marker.location <- data.obj$marker.location[-x.locale]
		}
		
	y.locale <- grep("Y", data.obj$chromosome, ignore.case = TRUE)
	if(length(y.locale) > 0){
		message("\nRemoving markers on the Y chromosome")
		gene <- gene[,-y.locale]
		data.obj$chromosome <- data.obj$chromosome[-y.locale]
		data.obj$marker.location <- data.obj$marker.location[-y.locale]
		}
	
	data.obj$geno <- gene
	
	n.phe <- dim(pheno)[2]
	n.gene <- dim(gene)[2]


	#first do the permutations to get the significance threshold
	#results will be compared to the significance threshold to 
	#determine
	if(verbose){
	cat("\nPerforming permutations to calculate significance threshold...\n")
	}
	
	data.obj <- genome.wide.threshold.1D(data.obj, n.perm = n.perm, alpha.for.pairs = alpha.for.pairs, alpha.for.covar = alpha.for.covar, scan.what = scan.what, verbose = verbose)
	
	#get the threshold for covariates
	covar.threshold <- data.obj$covar.thresh


	#if there are covariates specified, pull these out.
	#covariates must be coded as markers and contained in the
	#genotype matrix
	#first check the genotype and phenotype matrices for the
	#covariates. If the covariates are in the phenotype matrix
	#move them to the genotoype matrix. If you can't find the 
	#covariates, stop and warn the user.

	
	if(!is.null(covar)){

		covar.loc <- get.col.num(pheno, covar, warn = FALSE)
		if(length(covar.loc) > 0){
			stop("Phenotypic covariates should be created before running singlescan(). See create.covar() for help.")
			}
		
		covar.loc <- get.col.num(gene, covar, warn = FALSE)
		
		if(length(covar.loc) == 0){
			stop("I couldn't find the specified covariates. Please check the spelling, and make sure that covariates have been moved to the genotype matrix before this step.")
			}
		
		}else{
			covar.loc <- NULL
			}


		#This function gets regression statistics with a
		#covariate table
		get.stats <- function(phenotype, genotype.loc, covar.loc){

			#figure out if we are testing one of the covariates
			which.covar <- grep(genotype.loc, covar.loc)

			#if we are testing the covariate, take it out of
			#the covar.loc vector before testing
			if(length(which.covar) > 0){
				covar.loc <- covar.loc[-which.covar]
				}
				
				model <- lm(phenotype~gene[,c(covar.loc, genotype.loc)])
				
		
				#take the last line of coefficients.
				model.coef <- summary(model)$coefficients
				slope <- model.coef[dim(model.coef)[1],1]
				se <- model.coef[dim(model.coef)[1],2]
				t.stat <- abs(model.coef[dim(model.coef)[1],3])
				p.val <- model.coef[dim(model.coef)[1],4]
							
			#put together all the statistics we want to keep
			#we keep the absolute value of the t statistic,
			#the p value, and the covariate flag

			table.row <- c(slope, se, t.stat, p.val)
			return(table.row)
			}
	
	

	
	results.list <- vector(mode = "list", length = n.phe)	
	names(results.list) <- colnames(pheno)
	#==========================================

	for (i in 1:n.phe){
		
		#take out the response variable
		y <- pheno[,i]
				
		#apply the modeling function to each marker column
		results.table <- t(apply(matrix(c(1:dim(gene)[2]), nrow = 1), 2, function(x) get.stats(y, x, covar.loc)))
		colnames(results.table) <- c("slope", "se", "t.stat", "p.val")
		rownames(results.table) <- colnames(gene)
		results.list[[i]] <- results.table
		}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
	
		data.obj$singlescan.results <- results.list
	
		#calculate the covariate flags based on the oneD scan effects
		if(auto.covar.selection){
			data.obj <- get.covar(data.obj)
			}else{
			data.obj <- get.covar(data.obj, covar.thresh = Inf)	
			}
		data.obj$covar.thresh <- covar.threshold
		
		
		return(data.obj)
	
	}
