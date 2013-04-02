writeVariantInfluences <-
function(data.obj, pval.thresh = 0.05, filename = "Variant.Influences.csv", delim = ","){
	
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


	sig.var <- which(as.numeric(var.influences[,"p.adjusted"]) <= pval.thresh)
	
	if(length(sig.var) > 0){
		var.table <- var.influences[sig.var,]
		}else{
			var.table <- NULL
			}


	
	pheno.table <- NULL
	for(ph in pheno.names){
		sig.pheno <- which(as.numeric(pheno.results[[ph]][,"adj.p"]) <= pval.thresh)
		if(length(sig.pheno) > 0){
			pheno.section <- matrix(pheno.results[[ph]][sig.pheno,], nrow = length(sig.pheno))
			pheno.section  <- cbind(rep(ph, length(sig.pheno)), pheno.section)
			pheno.section[,1:2] <- pheno.section[,2:1] #awkward, but flip the source-target columns here.
			pheno.table <- rbind(pheno.table, pheno.section)
			}
		}
	
	
	#take out the raw t statistic column
	if(!is.null(pheno.table)){
		pheno.table <- pheno.table[,-5]
		colnames(pheno.table) <- colnames(var.table)
		}

	colnames(pheno.table) <- c("Source", "Target", "Effect", "SE", "|Effect|/SE", "P_empirical", "p.adjusted")

	final.table <- rbind(var.table, pheno.table)
	
	if(is.null(final.table)){
		final.table <- "No significant influences."
		}else{
			final.table <- final.table[order(final.table[,"|Effect|/SE"], decreasing = TRUE),]
			}

	write.table(final.table, file = filename, quote = FALSE, sep = delim, row.names = FALSE)
	
	invisible(final.table)	
}
