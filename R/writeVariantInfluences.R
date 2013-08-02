writeVariantInfluences <-
function(data.obj, p.or.q = 0.05, filename = "Variant.Influences.csv", delim = ",", mark.covar = FALSE, write.file = TRUE){
	
	var.influences <- data.obj$var.to.var.p.val
	pheno.results <- data.obj$max.var.to.pheno.influence
	
	if(is.null(var.influences)){
		stop("calc.p() must be run to calculate variant-to-variant influences.")
		}


	if(is.null(pheno.results)){
		stop("direct.influence() must be run to calculate variant-to-trait influences.")
		}

	
	if(data.obj$transform.to.phenospace){
		pheno.names <- colnames(data.obj$pheno)
		}else{
			pheno.names <- names(data.obj$pairscan.results)
			}
	num.pheno <- length(pheno.names)
	
	


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
			
	#replace 0 adjusted p values with a < string
	adj.p.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(final.table))))
	smallest.adj.p <- min(as.numeric(final.table[which(as.numeric(final.table[,adj.p.col]) > 0),adj.p.col]))
	zero.locale <- which(final.table[,adj.p.col] == 0)
	final.table[zero.locale,adj.p.col] <- paste("<", signif(smallest.adj.p, 3))

	if(mark.covar){
		covar.flags <- data.obj$covar.for.pairscan
		is.covar <- rownames(covar.flags)[unique(as.vector(unlist(apply(covar.flags, 2, function(x) which(x == 1)))))]
		covar.names <- data.obj$marker.names[match(is.covar, colnames(data.obj$geno))]
		covar.source.locale <- which(final.table[,1] %in% covar.names)
		covar.target.locale <- which(final.table[,2] %in% covar.names)
		if(length(covar.source.locale) > 0){
			final.table[covar.source.locale,1] <- paste(final.table[covar.source.locale,1], "*", sep = "")
			}
		if(length(covar.target.locale) > 0){
			final.table[covar.target.locale,2] <- paste(final.table[covar.target.locale,2], "*", sep = "")
			}
		}
	
	if(write.file){
		write.table(final.table, file = filename, quote = FALSE, sep = delim, row.names = FALSE)	
		invisible(final.table)
		}else{
			return(final.table)
			}
}
