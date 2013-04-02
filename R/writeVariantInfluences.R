writeVariantInfluences <-
function(data.obj, p.or.q = 0.05, filename = "Variant.Influences.csv", delim = ","){
	
	var.influences <- data.obj$var.to.var.p.val
	pheno.results <- data.obj$max.var.to.pheno.influence
	
	if(data.obj$transform.to.phenospace){
		pheno.names <- colnames(data.obj$pheno)
		}else{
			pheno.names <- names(data.obj$pairscan.results)
			}
	num.pheno <- length(pheno.names)
	
	
	if(is.null(var.influences)){
		stop("calc.p() must be run to calculate variant-to-variant influences.")
		}

	if(is.null(pheno.results)){
		stop("direct.influence() must be run to calculate variant-to-trait influences.")
		}


	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(var.influences))))

	sig.var <- which(as.numeric(var.influences[, var.sig.col]) <= p.or.q)
	
	
	if(length(sig.var) > 0){
		var.table <- var.influences[sig.var,,drop=FALSE]
		}else{
			var.table <- NULL
			}


	
	pheno.table <- NULL
	pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.results[[1]]))))
	for(ph in pheno.names){
		sig.pheno <- which(as.numeric(pheno.results[[ph]][,pheno.sig.col]) <= p.or.q)
		if(length(sig.pheno) > 0){
			pheno.section <- matrix(pheno.results[[ph]][sig.pheno,], nrow = length(sig.pheno))
			pheno.section  <- cbind(rep(ph, length(sig.pheno)), pheno.section)
			pheno.section[,1:2] <- pheno.section[,2:1] #awkward, but flip the source-target columns here.
			pheno.table <- rbind(pheno.table, pheno.section)
			}
		}
	
	if(length(pheno.table) > 0){
		pheno.table <- matrix(pheno.table, ncol = 8)
		#take out the raw t statistic column
		if(!is.null(pheno.table)){
			pheno.table <- pheno.table[,-5,drop=FALSE]
			colnames(pheno.table) <- colnames(var.table)
			}

		colnames(pheno.table) <- c("Source", "Target", "Effect", "SE", "|Effect|/SE", "P_empirical", colnames(pheno.results[[1]])[pheno.sig.col])
		}
		
	for(j in 1:2){
		marker.locale <- match(var.table[,j], colnames(data.obj$geno))
		var.table[,j] <- data.obj$marker.names[marker.locale]
		}

	marker.locale <- match(pheno.table[,1], colnames(data.obj$geno))
	pheno.table[,1] <- data.obj$marker.names[marker.locale]

	final.table <- rbind(var.table, pheno.table)
	
	if(is.null(final.table)){
		final.table <- "No significant influences."
		}else{
			final.table <- final.table[order(final.table[,"|Effect|/SE"], decreasing = TRUE),]
			}

	write.table(final.table, file = filename, quote = FALSE, sep = delim, row.names = FALSE)
	
	invisible(final.table)	
}
