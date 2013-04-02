plotSinglescan <-
function(data.obj, chr = NULL, traits = NULL, standardized = TRUE, mark.covar = TRUE, mark.chr = TRUE, plot.type = "h", overlay = FALSE, trait.colors = NULL, show.rejected.markers = FALSE, show.selected.markers = FALSE){
	
	if(show.rejected.markers && show.selected.markers){
		stop("show.rejected.markers and show.rejected.markers cannot both be TRUE.")
		}
	
	D1.results <- data.obj$singlescan.results
	marker.names <- data.obj$marker.names
	ind.markers <- data.obj$geno.for.pairscan
		
	if(is.null(D1.results)){
		stop("singlescan() must be run before plotting the results")
		}


	if(show.rejected.markers || show.selected.markers){
		if(is.null(ind.markers)){
			stop("select.markers.for.pairscan() must be run before plotting the selected markers.")
			}
		}

	if(is.null(chr)){
		chr <- unique(data.obj$chromosome)
		}

	if(is.null(traits)){
		traits <- names(D1.results)
		}
		

	covar.flags <- data.obj$covar.flags
	col.mat <- matrix(NA, nrow = dim(covar.flags)[1], ncol = dim(covar.flags)[2])

	if(!overlay){
		col.mat[covar.flags == 0] <- "black"
		if(mark.covar){
			col.mat[covar.flags == 1] <- "red"
			}else{
			col.mat[covar.flags == 1] <- "black"	
			}
		}else{
		if(is.null(trait.colors)){
			trait.colors <- c("black", "blue", "purple", "darkgreen")
			}
		if(length(trait.colors) < length(traits)){
		 	trait.colors <- rep(trait.colors, length(traits)/4)
		 	trait.colors <- trait.colors[1:length(traits)]
			}
		for(i in 1:length(traits)){
			col.mat[covar.flags[,i] == 0, i] <- trait.colors[i]
			if(mark.covar){
				col.mat[covar.flags[,i] == 1,i] <- "red"
				}else{
				col.mat[covar.flags[,i] == 1,i] <- trait.colors[i]
				}
			}
			
		}


	results.rows <- which(data.obj$chromosome %in% chr)
	results.el <- which(names(D1.results) %in% traits)
	results.to.plot <- NULL
	for(r in results.el){
		if(standardized){
			results.to.plot <- cbind(results.to.plot, D1.results[[r]][results.rows,"t.stat"])
			}else{
			results.to.plot <- cbind(results.to.plot, D1.results[[r]][results.rows,"slope"])
			}
		}


	final.cols <- col.mat[results.rows, results.el]

	colnames(results.to.plot) <- names(D1.results)[results.el]


		if(show.rejected.markers || show.selected.markers){
			if(show.selected.markers){
				ind.locale <- which(rownames(results.to.plot) %in% colnames(ind.markers) )
				}else{
				ind.locale <- which(!rownames(results.to.plot) %in% colnames(ind.markers))
				}
			}


	pairscan.threshold <- data.obj$pairscan.thresh
	covar.threshold <- data.obj$covar.thresh

	if(!overlay){
		layout.mat <- get.layout.mat(length(results.el), "upright")
		}else{
		layout.mat <- matrix(1, 1, 1)
		}

	layout(layout.mat)	
	for(p in 1:length(results.to.plot[1,])){
		# dev.new(width = plot.width, height = plot.height)
		pheno.res <- results.to.plot[,p]

		if(standardized){
			if(!overlay){
				all.vals <- c(pheno.res, pairscan.threshold, covar.threshold, 0)		
				}else{
				all.vals <- c(results.to.plot, pairscan.threshold, covar.threshold, 0)	
					}
			}else{
				if(!overlay){
					all.vals <- pheno.res
					}else{
					all.vals <- results.to.plot	
					}
			}

		#create the window
		if(p == 1 || !overlay){
			par(mar = c(3, 4, 3, 2) + 0.1)
			plot.new()
			

			plot.window(xlim = c(0, length(pheno.res)), ylim = c(min(all.vals), max(all.vals[is.finite(all.vals)])))
			
		
			#shade the chromosome regions
			if(mark.chr){
				markers.used.locale <- which(colnames(data.obj$geno) %in% rownames(results.to.plot))
				chr.id <- data.obj$chromosome[markers.used.locale]
				par(xpd = TRUE)
				for(ch in 1:length(chr)){
					x.min <- min(which(chr.id == chr[ch])); x.max <- max(which(chr.id == chr[ch]))
					if(ch %% 2 == 1){
						polygon(x = c(x.min, x.min, x.max, x.max), y = c(min(all.vals), max(all.vals), max(all.vals), min(all.vals)), col = "lightgray", border = NA)
						}
					if(chr[ch] == 0){
						text(x = x.max, y = min(all.vals)-((max(all.vals)-min(all.vals))*0.05), labels = "Cov.", cex = 0.5, adj = 0)
						}else{
						text(x = mean(c(x.min, x.max)), y = min(all.vals)-((max(all.vals)-min(all.vals))*0.05), labels = chr[ch], cex = 0.5)
						}
					}
				par(xpd = FALSE)
				}
		
				# axis(1, at = 1:length(pheno.res), labels = FALSE)
				abline(h = 0)
			    # lbl <- marker.names


			if(show.selected.markers){
				ind.locale <- which(rownames(results.to.plot) %in% colnames(ind.markers) )
				}else{
				ind.locale <- which(!rownames(results.to.plot) %in% colnames(ind.markers))
				}


			if(standardized){

				#add the marker labels to the x axis	    
			    # text(1:length(pheno.res), par("usr")[3] - 0.25, srt = 90, adj = 1, labels = NULL, xpd = TRUE, cex = 0.5)
				#add lines to indicate the significance thresholds
				abline(h = pairscan.threshold, lty = 1)
				abline(h = covar.threshold, lty = 2)
				par(xpd = TRUE)
				text(x = length(marker.names)*1.05, y = pairscan.threshold, labels = paste("p =", data.obj$alpha.for.pairs), cex = 0.5, adj = 0)
				text(x = length(marker.names)*1.05, y = covar.threshold, labels = paste("p =", data.obj$alpha.for.covar), cex = 0.5, adj = 0)
				par(xpd = FALSE)
								
				if(!overlay){
					mtext(colnames(results.to.plot)[p], cex = 2)
					mtext(paste(colnames(results.to.plot)[p], "[|Eff|/se]", sep = " "), side = 2, line = 2.5)
					}else{
					mtext("|Eff|/se", side = 2, line = 2.5)	
					}
				axis(2)

				}else{
				
				abline(h = 0)
				if(!overlay){
					mtext(colnames(results.to.plot)[p], cex = 2)
					mtext(paste(colnames(results.to.plot)[p], "[Eff]", sep = " "), side = 2, line = 2.5)	
					}else{
					mtext("Eff", side = 2, line = 2.5)		
					}
				
				axis(2)

				# text(1:length(pheno.res), par("usr")[3] - abs((par("usr")[3])*0.1), srt = 90, adj = 1, xpd = TRUE, cex = 0.5)
				}
				
			if(show.selected.markers || show.rejected.markers){
				points(ind.locale, (pheno.res[ind.locale]+(max(pheno.res)*0.02)), col = "red", pch = "*")
				}

		
			if(p == 1){
				par(xpd = TRUE)
				if(mark.covar){
					if(plot.type == "p" || plot.type == "b"){
						legend(x = (0-(length(marker.names)*0.04)), y = max(all.vals)*1.15, legend = "covariate", pch = 16, col = "red", cex = 0.7)
						}
					if(plot.type == "h"){
						legend(x = (0-(length(marker.names)*0.04)), y = max(all.vals)*1.15, legend = "covariate", lty = 1, col = "red", cex = 0.7)
						}
					}
				if(show.selected.markers){
					legend(x = (length(marker.names)*0.94), y = max(all.vals)*1.15, legend = "selected", pch = "*", col = "red", cex = 0.7)
					}
				if(show.rejected.markers){
					legend(x = (length(marker.names)*0.94), y = max(all.vals)*1.15, legend = "rejected", pch = "*", col = "red", cex = 0.7)
					}
				par(xpd = FALSE)
				}
			}
			
		if(standardized){
			#plot the effect sizes
			points(pheno.res, type = plot.type, col = col.mat[,p], pch = 16)
			if(!overlay){
				mtext(colnames(results.to.plot)[p], cex = 2)
				mtext(paste(colnames(results.to.plot)[p], "[|Eff|/se]", sep = " "), side = 2, line = 2.5)
				axis(2)
				}
			}else{
			points(pheno.res, type = plot.type, col = col.mat[,p], pch = 16)
			if(!overlay){
				abline(h = 0)
				mtext(colnames(results.to.plot)[p], cex = 2)
				mtext(paste(colnames(results.to.plot)[p], "[Eff]", sep = " "), side = 2, line = 2.5)
				axis(2)
				}
			}
			
		}
		
		if(overlay){
			par(xpd = TRUE)
			if(show.rejected.markers || show.selected.markers){
				legend(x = (length(marker.names)*0.82), y = max(all.vals)*1.15, legend = traits, pch = 16, col = trait.colors[1:length(traits)], cex = 0.7)
				}else{
				legend(x = (length(marker.names)*0.94), y = max(all.vals)*1.15, legend = traits, pch = 16, col = trait.colors[1:length(traits)], cex = 0.7)	
				}
			par(xpd = FALSE)
		}
		
			

}
