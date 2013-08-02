one.pairscan <-
function(phenotype.vector, genotype.matrix, covar.vector, pairs.matrix, n.perm = 0, verbose = FALSE){
			
		# require(Matrix)
		
		#============================================================================
		# check to see that the covariates are not redundant
		#============================================================================
		use.covars <- which(covar.vector == 1)
		if(length(use.covars) > 0){
			cov.mat <- matrix(genotype.matrix[,as.numeric(use.covars)], ncol = length(use.covars))
			cor.mat <- cor(cov.mat)
			diag(cor.mat) <- 0
			perfect.cor <- which(abs(signif(cor.mat, 2)) == 1)
			if(length(perfect.cor) > 0){
				stop("Check the covariates. There appears to be at least one pair of redundant covariates.")
				}
			}
		#============================================================================
		
			
		#============================================================================
		#internal functions
		#============================================================================
		
			
		get.model.results <- function(marker.names, m1, m2){
				
				model.effects <- matrix(0, ncol = (4+length(covar.vector)), nrow = 1)	#initialize a row to hold the model statistics for the pair. 
															 							#It has enough colums for the intercept (1), max number
															 							#of covariates (length(covar.vector)),
															 							#individual markers (2), and the interaction between 
															 							#the markers (1)
															 							#It starts out filled with 0s because we fix all beta
															 							#coefficients for non-covariates at 0
				model.se <- matrix(0, ncol = (4+length(covar.vector)), nrow = 1)
				model.cov.results <- matrix(NA, ncol = 9, nrow = 1) #a table to hold a linearized covariance matrix of the model for each pair
				# model.cov.results <- matrix(NA, ncol = 6, nrow = 1)
				
				#check to see if any of the markers we are testing
				#in the pair are listed as covariates. If they are,
				#we exclude them from the covariate list.
				#use.covars ends up as a vector of marker columns
				#to use as covariates.
				use.covars <- which(covar.vector == 1)
				for(i in 1:length(marker.names)){
					marker.locale <- which(names(use.covars) == marker.names[i])
					if(length(marker.locale) > 0){
						use.covars <- use.covars[-marker.locale]
					}
				}
				
				#create a vector that will dictate where we put
				#the resulting coefficients from the linear model
				#there will be one spot for the intercept, one
				#each for each covariate, and one each for marker1
				#marker2 and marker1:marker2
				coefficient.locale <- as.vector(c(1, (use.covars+1), (dim(model.effects)[2]-2):dim(model.effects)[2]))
			
				#pull out the genotype data for each marker for ease of 
				#reading the code later
			
				design.mat <- cbind(rep(1, length(m1)), genotype.matrix[,use.covars], m1, m2)
				missing.rows <- which(is.na(rowSums(design.mat)))
				if(length(missing.rows) > 0){
					design.mat <- design.mat[-missing.rows,]
					}
				rank.check <- rankMatrix(design.mat)[[1]]

				if(rank.check < dim(design.mat)[2]){
					return(NULL)
					}

				#Do the linear regression with the covariates, the two markers 
				#individually, and the interaction between the two markers
				#put the covars first so the effects of the markers are calculated
				#after the covariates have been taken into account
				if(length(use.covars) > 0){
					model <- lm(phenotype.vector ~ genotype.matrix[,use.covars] + m1 + m2 + m1:m2, model = FALSE, x = FALSE, y = FALSE, qr = TRUE)
					}else{
						model <- lm(phenotype.vector ~ m1 + m2 + m1:m2, model = FALSE, x = FALSE, y = FALSE, qr = TRUE)	
						}
				model.summ <- summary(model)
				
				#occasionally permutations result in non linearly dependent matrices.
				#if this is the case, return NULL. This triggers the permutation
				#script to generate another permutation.
				if(length(which(is.na(coefficients(model)))) > 0){ 
					return(NULL)
					}
			
				#we need to put the entries in the correct spots in the coefficient matrix
				#cat("Phenotype = ",p,"Marker Pair = ",m,"\n")
				
				model.effects[,coefficient.locale] <- coef(model)
				
				model.se[,coefficient.locale] <- model.summ$coefficients[,"Std. Error"]

				#calculate the covariance matrix of the model parameters
				model.cov.mat <- vcov(model)
				dim.mat <- dim(model.cov.mat)[1] #get the dimensions of the matrix. We want the last three rows and last three columns (the covariance matrix for m1, m2, and m1:m2)
				cov.mat <- model.cov.mat[(dim.mat-2):dim.mat, (dim.mat-2):dim.mat]
				model.cov.results <- as.vector(cov.mat)

				results <- list(model.effects, model.se, model.cov.results)
				names(results) <- c("model.effects", "model.se", "model.cov")
				return(results)
				}
		#============================================================================
		#end of internal functions
		#============================================================================
		
		# each.iter.time <- rep(NA, n.pairs)
		#we can calculate the number of genes
		#from the input data
		n.gene <- dim(genotype.matrix)[2]
			
		n.pairs <- dim(pairs.matrix)[1]
	
		#the locations of the 1's indicate which markers were
		#significantly associated with the phenotype
		#these markers will be put into the model
		#as covariates

		all.model.effects <- matrix(NA, nrow = n.pairs, ncol = (1+length(covar.vector)+2+1))
									#initialize a matrix to hold the model statistics for the pair. 
		 							#It has enough columns for the intercept (1), max number
		 							#of covariates (length(covar.vector)),
		 							#individual markers (2), and the interaction between 
		 							#the markers (1)
		 							#It starts out filled with 0s because we fix all beta
		 							#coefficients for non-covariates at 0

		all.model.se <- matrix(NA, nrow = n.pairs, ncol = (1+length(covar.vector)+2+1))

		all.model.cov <- matrix(NA, nrow = n.pairs, ncol = 9) #a table to hold a linearized covariance matrix of the model for each pair

		#also make variables to hold the permutation results
		all.model.effects.perm <- matrix(NA, nrow = n.pairs*n.perm, ncol = (1+length(covar.vector)+2+1))
		all.model.se.perm <- matrix(NA, nrow = n.pairs*n.perm, ncol = (1+length(covar.vector)+2+1))

		all.model.cov.perm <- matrix(NA, nrow = n.pairs*n.perm, ncol = 9)

		marker.pairs.used.perm <- matrix(NA, nrow = n.pairs*n.perm, ncol = 2)
		
		perm.position <- 1

		#for each of the marker pair combinations
		for (m in (1:length(pairs.matrix[,1]))){
			# begin.time <- Sys.time()
			# cat(m, "\n")
			if(verbose){
				report.progress(m, n.pairs, percent.text = 10, percent.dot = 2)
				}
				
			#get the marker identities
			markers <- as.vector(pairs.matrix[m,])
			marker1 <- genotype.matrix[,markers[1]]
			marker2 <- genotype.matrix[,markers[2]]

			marker.pair.results <- get.model.results(marker.names = markers, m1 = marker1, m2 = marker2)
			
			if(!is.null(marker.pair.results)){
				all.model.effects[m,] <- marker.pair.results[[1]]
				all.model.se[m,] <- marker.pair.results[[2]]
				all.model.cov[m,] <- marker.pair.results[[3]]
				
				#run ther permutations if specified
				if(n.perm > 0){
					p <- 0
					while(p < n.perm){
						# cat("\tDoing permutations\n")						
						#shuffle the markers in tandem and run them through the test again
						rnd.order <- sample(1:length(marker1), length(marker1))
						marker1.shuff <- marker1[rnd.order]
						marker2.shuff <- marker2[rnd.order]
						
						marker.shuff.results <- get.model.results(markers, marker1.shuff, marker2.shuff)
						#only record data and increment the permutation count if we got a good permutation
						if(!is.null(marker.shuff.results)){
							marker.pairs.used.perm[perm.position,] <- c(markers[1], markers[2])
							all.model.effects.perm[perm.position,] <- marker.shuff.results[[1]]
							all.model.se.perm[perm.position,] <- marker.shuff.results[[2]]
							all.model.cov.perm[perm.position,] <- marker.shuff.results[[3]]
							p <- p + 1
							perm.position <- perm.position + 1
							}

						} #end looping through permutations
					
					} #end case for if permutations are asked for
				}#end case for marker.pair.results being NULL
				# end.time <- Sys.time()
				# each.iter.time[m] <- as.numeric(end.time-begin.time)
	    	} #end looping through marker pairs
	      
	      
	      	#assign column names to the results tables
			#the column names represent the possible beta
			#coefficients we can get. The intercept, all
			#possible covariates, marker1 and marker2,
			#and the interaction marker1:marker2
			column.names <- c("(Intercept)", colnames(genotype.matrix), "marker1", "marker2", "marker1:marker2")
			colnames(all.model.effects) <- colnames(all.model.se) <- column.names

		     untested.pairs <- which(is.na(rowSums(all.model.effects)))
		     if(length(untested.pairs) > 0){
		     	all.model.effects <- all.model.effects[-untested.pairs,]
		     	all.model.se <- all.model.se[-untested.pairs,]
		     	all.model.cov <- all.model.cov[-untested.pairs,]
		     	pairs.matrix <- pairs.matrix[-untested.pairs,]
		     	}
		      
		      
			if(!is.null(all.model.effects.perm)){
				colnames(all.model.effects.perm) <- colnames(all.model.se.perm) <- column.names
				untested.pairs.perm <- which(is.na(rowSums(all.model.effects.perm)))
				if(length(untested.pairs.perm) > 0){
			     	all.model.effects.perm <- all.model.effects.perm[-untested.pairs.perm,]
			     	all.model.se.perm <- all.model.se.perm[-untested.pairs.perm,]
			     	all.model.cov.perm <- all.model.cov.perm[-untested.pairs.perm,]
		    	 	marker.pairs.used.perm <- marker.pairs.used.perm[-untested.pairs.perm,]
					}	
				}
	      
			#add the marker pair names to the results tables and name the columns
	      	colnames(pairs.matrix) <- c("marker1", "marker2")
	      	pairs.matrix <- apply(pairs.matrix, 2, as.numeric)
			final.effects.table <- cbind(pairs.matrix, all.model.effects)
			final.se.table <- cbind(pairs.matrix, all.model.se)
			
			if(!is.null(all.model.effects.perm)){
				marker.pairs.used.perm <- apply(marker.pairs.used.perm, 2, as.numeric)
		      	colnames(marker.pairs.used.perm) <- c("marker1", "marker2")
				final.effects.table.perm <- cbind(marker.pairs.used.perm, all.model.effects.perm)
				final.se.table.perm <- cbind(marker.pairs.used.perm, all.model.se.perm)	
				}

			final.cov.table <- all.model.cov
			
			if(!is.null(all.model.effects.perm)){
				final.cov.table.perm <- all.model.cov.perm
				}

			phenotype.results <- list(final.effects.table, final.se.table, final.cov.table)
			names(phenotype.results) <- c("pairscan.effects", "pairscan.se", "model.covariance")
	      	
	      	if(!is.null(all.model.effects.perm)){
				phenotype.perm.results <- list(final.effects.table.perm, final.se.table.perm, final.cov.table.perm)
				names(phenotype.perm.results) <- c("pairscan.effects.perm", "pairscan.se.perm", "model.covariance.perm")
	      		}else{
	      			phenotype.perm.results <- NULL
	      			}
	      	
	      	
	      	final.results <- list(phenotype.results, phenotype.perm.results)
	      	names(final.results) <- c("pairscan.results", "pairscan.perm")

			if(verbose){
				cat("\n") #make sure the prompt is on the next line at the end of everything
				}


	      	return(final.results)


	}
