plotSinglescan <-
function(data.obj, chr = NULL, traits = NULL, show.alpha.values = NULL, standardized = TRUE, show.marker.labels = FALSE, mark.covar = TRUE, mark.chr = TRUE, plot.type = "h", overlay = FALSE, trait.colors = NULL, show.rejected.markers = FALSE, show.selected.markers = FALSE){
	
	
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


		all.alpha <- data.obj$alpha
		
		if(!is.null(show.alpha.values)){
			alpha.to.include <- which(all.alpha %in% sort(show.alpha.values))
			if(length(alpha.to.include) < length(show.alpha.values)){
				cant.find <- setdiff(show.alpha.values, all.alpha)
				warning("The following alpha values were not calculated by singlescan():\n")
				cat(cant.find, sep = "\n")
				}
			}else{
			alpha.to.include <- 1:length(data.obj$alpha)
			}
			
		alpha.thresholds <- rep(NA, length(alpha.to.include))
		for(a in 1:length(alpha.to.include)){
			alpha.thresholds[a] <- data.obj$alpha.thresh[[alpha.to.include[a]]]
			}
		
		

	if(!overlay){
		layout.mat <- get.layout.mat(length(results.el), "upright")
		}else{
		layout.mat <- matrix(1, 1, 1)
		}

	layout(layout.mat)	
	for(p in 1:length(results.to.plot[1,])){

		#figure out the axix label
		if(!overlay){
			if(standardized){
				y.label <- paste(colnames(results.to.plot)[p], "[|Eff|/se]", sep = " ")
				}else{
				y.label <- paste(colnames(results.to.plot)[p], "Eff", sep = " ")	
				}
			}else{
			if(standardized){
				y.label <- "[|Eff|/se]"
				}else{
				y.label <- "Eff"
				}				
			}

		pheno.res <- results.to.plot[,p]

		if(standardized){
			if(!overlay){
				all.vals <- c(pheno.res, unlist(alpha.thresholds), 0)		
				}else{
				all.vals <- c(results.to.plot, unlist(alpha.thresholds), 0)
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
			# if(!show.marker.labels){
				markers.used.locale <- which(colnames(data.obj$geno) %in% rownames(results.to.plot))
				chr.id <- data.obj$chromosome[markers.used.locale]
				par(xpd = TRUE)
				for(ch in 1:length(chr)){
					x.min <- min(which(chr.id == chr[ch])); x.max <- max(which(chr.id == chr[ch]))
					if(ch %% 2 == 1 && mark.chr){ #shade chromosome regions if mark.chr is TRUE
						polygon(x = c(x.min, x.min, x.max, x.max), y = c(min(all.vals), max(all.vals), max(all.vals), min(all.vals)), col = "lightgray", border = NA)
						}

					if(!show.marker.labels){ #plot the chromosome labels if we are not plotting marker labels
						if(chr[ch] == 0){
							text(x = x.max, y = min(all.vals)-((max(all.vals)-min(all.vals))*0.05), labels = "Cov.", cex = 0.5, adj = 0, font = 2)
							}else{
							text(x = mean(c(x.min, x.max)), y = min(all.vals)-((max(all.vals)-min(all.vals))*0.05), labels = chr[ch], cex = 0.5, font = 2)
							}
						}
					}
				par(xpd = FALSE)
				# }
		
				abline(h = 0)
				
				if(show.marker.labels){
					if(standardized){
						if(mark.chr){y.val <- max(all.vals)*-0.08}else{y.val <- max(all.vals)*-0.05}
						}else{
						if(mark.chr){y.val <- min(all.vals) - (max(all.vals)*0.08)}else{y.val <- min(all.vals) - (max(all.vals)*0.05)}
						
						}
					# axis(1, at = 1:length(pheno.res), labels = FALSE)			
				    par(xpd = TRUE)
					text(x = 1:length(pheno.res), y = y.val, labels = marker.names, srt = 90, cex = 0.5, adj = 1)
					par(xpd = FALSE)
					}

			if(show.selected.markers){
				ind.locale <- which(rownames(results.to.plot) %in% colnames(ind.markers) )
				}else{
				ind.locale <- which(!rownames(results.to.plot) %in% colnames(ind.markers))
				}


			if(standardized){

				#add the lines for the significance thresholds
				for(a in 1:length(alpha.thresholds))				
					points(x = 1:length(marker.names), y = rep(alpha.thresholds[[a]], length(marker.names)), type = "l", lty = a)
					}

				#add labels for significance lines
				par(xpd = TRUE)
				for(a in 1:length(data.obj$alpha)){
					text(x = length(marker.names)*1.02, y = alpha.thresholds[[a]], labels = paste("p =", data.obj$alpha[a]), cex = 0.5, adj = 0)
					}
									
				par(xpd = FALSE)
								
				if(!overlay){
					mtext(colnames(results.to.plot)[p], cex = 2)
					mtext(y.label, side = 2, line = 2.5)
					}else{
					mtext(y.label, side = 2, line = 2.5)	
					}
		
				axis(2)
				abline(h = 0)
				
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
			
			
		marker.chr <- data.obj$chromosome[which(names(pheno.res) %in% colnames(data.obj$geno))]
		if(standardized){
			#plot the effect sizes
			for(ch in chr){
				chr.locale <- which(marker.chr == ch)
				if(ch != 0){
					points(x = chr.locale, y = pheno.res[chr.locale], type = plot.type, col = col.mat[chr.locale,p], pch = 16)
					}else{
					points(x = chr.locale, y = pheno.res[chr.locale], type = "h", col = col.mat[chr.locale,p], pch = 16)	
					}
				}
			
			}else{
			
			for(ch in chr){
				chr.locale <- which(marker.chr == ch)
				if(ch != 0){
					points(x = chr.locale, y = pheno.res[chr.locale], type = plot.type, col = col.mat[,p], pch = 16)
					}else{
					points(x = chr.locale, y = pheno.res[chr.locale], type = "h", col = col.mat[,p], pch = 16)	
					}
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
