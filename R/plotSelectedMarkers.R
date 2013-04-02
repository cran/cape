plotSelectedMarkers <-
function(data.obj, mark.which = c("selected", "rejected"), chr = NULL, traits = NULL, standardized = TRUE, mark.covar = FALSE, mark.chr = FALSE){
	
	D1.results <- data.obj$singlescan.results
	ind.markers <- data.obj$geno.for.pairscan
	
	if(length(grep("rej", mark.which)) > 0){
		mark.which <- "rejected"
		}else{
		mark.which <- "selected"
		}
	
	if(is.null(D1.results)){
		stop("singlescan() must be run before plotting the selected markers.")
		}

	if(is.null(ind.markers)){
		stop("select.markers.for.pairscan() must be run before plotting the selected markers.")
		}

	if(is.null(chr)){
		chr <- unique(data.obj$chromosome)
		}

	if(is.null(traits)){
		traits <- names(D1.results)
		}
	

	covar.flags <- data.obj$covar.flags
	col.mat <- matrix(NA, nrow = dim(covar.flags)[1], ncol = dim(covar.flags)[2])
	col.mat[covar.flags == 0] <- "black"
	if(mark.covar){
		col.mat[covar.flags == 1] <- "red"
		}else{
		col.mat[covar.flags == 1] <- "black"	
		}
			

	#Get the marker we'll be plotting
	results.rows <- which(data.obj$chromosome %in% chr)
	found.chr <- which(chr %in% data.obj$chromosome)
	if(length(found.chr) < length(chr)){
		if(length(found.chr) > 0){
			not.found <- chr[-found.chr]
			}else{
			not.found <- chr
			}
		message("\nI couldn't find the following chromosomes:")
		cat(not.found, sep = "\n")
		return()
		}
	
	
	results.el <- which(names(D1.results) %in% traits)
	
	if(length(results.el) < length(traits)){
		if(length(results.el) > 0){
			not.found <- traits[-results.el]
			}else{
			not.found <- traits
			}
		message("I couldn't find the following traits:")
		cat(not.found, sep = "\n")
		return()
		}
		
	
	results.to.plot <- NULL
	for(r in results.el){
		if(standardized){
			results.to.plot <- cbind(results.to.plot, D1.results[[r]][results.rows,"t.stat"])
			}else{
			results.to.plot <- cbind(results.to.plot, D1.results[[r]][results.rows,"slope"])
			}
		}
		
	final.cols <- matrix(col.mat[results.rows, results.el], ncol = length(results.el))

	if(mark.which == "selected"){
		ind.locale <- which(rownames(results.to.plot) %in% colnames(ind.markers) )
		}else{
		ind.locale <- which(!rownames(results.to.plot) %in% colnames(ind.markers))
		}

	colnames(results.to.plot) <- names(D1.results)[results.el]

	pairscan.threshold <- data.obj$pairscan.thresh
	covar.threshold <- data.obj$covar.thresh

		
	layout.mat <- get.layout.mat(length(results.el), "upright")
	
	layout(layout.mat)
	for(p in 1:length(results.to.plot[1,])){
		# dev.new(width = plot.width, height = plot.height)
		pheno.res <- results.to.plot[,p, drop = FALSE]

		if(standardized){
			all.vals <- c(pheno.res, pairscan.threshold, covar.threshold, 0)		
			}else{
			all.vals <- pheno.res
			}

		#create the window
		par(mar = c(3, 4, 3, 2) + 0.1)
		plot.new()
		plot.window(xlim = c(0, length(pheno.res)), ylim = c(min(all.vals), max(all.vals)))

		#shade the chromosome regions
		if(mark.chr){
			markers.used.locale <- which(colnames(data.obj$geno) %in% rownames(results.to.plot))
			chr.id <- data.obj$chromosome[markers.used.locale]
			for(ch in seq(1,length(chr), 2)){
				x.min <- min(which(chr.id == chr[ch])); x.max <- max(which(chr.id == chr[ch]))
				polygon(x = c(x.min, x.min, x.max, x.max), y = c(min(all.vals), max(all.vals), max(all.vals), min(all.vals)), col = "lightgray", border = NA)
				}
			}

		if(standardized){
			#plot the effect sizes
			points(pheno.res, type="h", col = col.mat[,p])
			mtext(colnames(results.to.plot)[p], cex = 2)
			mtext(paste(colnames(results.to.plot)[p], "[|Eff|/se]", sep = " "), side = 2, line = 2.5)
			axis(2)
			}else{
			points(pheno.res, type="h", col = col.mat[,p])		
			abline(h = 0)
			mtext(colnames(results.to.plot)[p], cex = 2)
			mtext(paste(colnames(results.to.plot)[p], "[Eff]", sep = " "), side = 2, line = 2.5)
			axis(2)
			}

					
		if(!is.null(mark.which)){
			points(ind.locale, (pheno.res[ind.locale]+(max(pheno.res)*0.02)), col = "red", pch = "*")
			}

		
		if(p == 1){
			if(!is.null(mark.which)){
				if(mark.which == "selected"){
					legend.label = "Selected Marker"
					}else{
					legend.label = "Discarded Marker"
					}
				legend("topleft", legend = legend.label, col = "red", pch = "*", cex = 0.7)
				}
			
				if(mark.covar){
					legend("topright", legend = "Covariate", lty = 1, col = "red", cex = 0.7)
					}

			}
		

		
		axis(1, at = 1:length(pheno.res), labels = FALSE)
	   
	    lbl <- rownames(results.to.plot)
	
		if(standardized){
			#add the marker labels to the x axis	    
		    text(1:length(pheno.res), par("usr")[3] - 0.25, srt = 90, adj = 1, labels = lbl, xpd = TRUE, cex = 0.5)
			#add lines to indicate the significance thresholds
			abline(h = pairscan.threshold, lty = 1)
			abline(h = covar.threshold, lty = 2)
			}
	
		}

	
	



}
