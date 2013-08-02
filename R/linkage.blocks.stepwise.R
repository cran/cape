linkage.blocks.stepwise <-
function(data.obj, collapse.linked.markers = TRUE, r.thresh = 0.7, plot.blocks = FALSE){
	
	net.data <- data.obj$var.to.var.p.val
	pheno.net.data <- data.obj$max.var.to.pheno.influence
	
	
	if(length(net.data) == 0){
		stop("calc.p() must be run to calculate variant-to-variant influences.")
		}
		
	
	# #get the markers with significant influences on other markers
	# #and on phenotypes
	# var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(net.data))))
	# sig.markers <- net.data[which(net.data[, var.sig.col] <= p.or.q), 1:2]
	
	# sig.pheno <- vector(mode = "list", length = length(pheno.net.data))
	# pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.net.data[[1]]))))
	# for(ph in 1:length(pheno.net.data)){
		# sig.pheno[[ph]] <- pheno.net.data[[ph]][which(pheno.net.data[[ph]][,pheno.sig.col] <= p.or.q),1]
		# }

	# if(length(c(sig.markers, unlist(sig.pheno))) == 0){
		# stop("There are no significant markers at p = ", p.or.q)
		# }

	# #find all the markers with significant p values
	# #and sort them
	# u_markers <- unique(c(sig.markers, unlist(sig.pheno)))
	# marker.locale <- match(u_markers, colnames(data.obj$geno))
	# sorted.markers <- u_markers[order(marker.locale)]
	# marker.chr <- data.obj$chromosome[sort(marker.locale)]
	# u_chr <- unique(marker.chr)

	#find all the chromosomes that were used in the pairwise scan and sort them
	used.markers <- colnames(data.obj$geno)[which(colnames(data.obj$geno) %in% colnames(data.obj$geno.for.pairscan))]
	marker.chr <- data.obj$chromosome[which(colnames(data.obj$geno) %in% colnames(data.obj$geno.for.pairscan))]
	u_chr <- sort(unique(as.numeric(marker.chr)))
	
	#========================================================================================
	# internal functions
	#========================================================================================
	

	#if we are not collapsing the markers into blocks, 
	#or a chromosome only has one marker, or we are
	#adding the covariate chromosome, just add all 
	#individual markers to the list of blocks.
	add.ind.markers <- function(link.blocks, ch, chr.markers){
		chr.block.num <- 1
		if(is.null(link.blocks[[1]])){
			num.blocks = 1
			}else{
			num.blocks = length(link.blocks) + 1	
			}
		for(i in 1:length(chr.markers)){
			link.blocks[[num.blocks]] <- chr.markers[i]
			if(ch == 0){
				names(link.blocks)[num.blocks] <- data.obj$marker.names[which(colnames(data.obj$geno) == chr.markers[i])]
				}else{
				names(link.blocks)[num.blocks] <- paste("Chr", ch, "_", chr.block.num, sep = "")
				}
			num.blocks <- num.blocks + 1
			chr.block.num <- chr.block.num + 1
			}
		return(link.blocks)
		}
	

	#get the recombination data for the markers on a given chromosome
	get.chr.cor <- function(chr.markers){		
		just.markers <- unique(chr.markers)
		marker.locale <- which(colnames(data.obj$geno) %in% just.markers)
		chr.geno <- data.obj$geno[,marker.locale]
		chr.cor <- cor(chr.geno, use = "complete.obs")
		chr.cor[lower.tri(chr.cor, diag = FALSE)] <- NA
		return(chr.cor)
		}
		

	#I'm trying by traveling across the recomb.rf row, and stopping
	#at the user-imposed rf.min
	search.blocks.forward <- function(all.cor, r.thresh){
		start.row <- 1
		num.blocks <- 1
		new.blocks <- vector(mode = "list", length = 1)
		while(start.row <= dim(all.cor)[1]){
			block.markers <- which(all.cor[start.row,] >= r.thresh)
			if(length(block.markers) == 0){
				new.blocks[[num.blocks]] <- rownames(all.cor)[start.row]
				start.row <- start.row + 1
				}else{
				block <- colnames(all.cor)[min(block.markers):max(block.markers)]
				new.blocks[[num.blocks]] <- block
				start.row <- max(block.markers) + 1
				}
			names(new.blocks)[num.blocks] <- paste("Chr", ch, "_", num.blocks, sep = "")
			num.blocks <- num.blocks + 1
			}
		return(new.blocks)
		}

	
	add.chr.blocks <- function(link.blocks, new.blocks){
		chr.block.num <- 1
		if(is.null(link.blocks[[1]])){
			num.blocks = 1
			}else{
			num.blocks = length(link.blocks) + 1	
			}		
		for(i in 1:length(new.blocks)){
			link.blocks[[num.blocks]] <- new.blocks[[i]]
			names(link.blocks)[num.blocks] <- names(new.blocks)[i]
			num.blocks <- num.blocks + 1
			}
		return(link.blocks)
		}
	
	
	#========================================================================================
	# end internal functions
	#========================================================================================


	if(plot.blocks){pdf(paste("Recomb.Images.r.thresh.", r.thresh, ".pdf", sep = ""))}
	#go through each chromosome separately and find the linkage blocks on each chromosome
	link.blocks <- vector(mode = "list", length = 1)
	num.blocks <- 1
	for(ch in u_chr){
		chr.markers <- used.markers[which(as.numeric(marker.chr) == as.numeric(ch))]
		
		if(!collapse.linked.markers || length(chr.markers) == 1 || ch == 0){
			link.blocks <- add.ind.markers(link.blocks, ch, chr.markers)
			}else{
			# all.recomb <- get.chr.recomb(chr.markers)
			all.cor <- get.chr.cor(chr.markers)
			forward.blocks <- search.blocks.forward(all.cor, r.thresh)
			# backward.blocks <- search.blocks.backward(all.cor, r.thresh) 
			link.blocks <- add.chr.blocks(link.blocks, forward.blocks)
			
			if(plot.blocks){
				# all.recomb[which(all.recomb >= rf.min)] <- NA
				image(x = 0:(dim(all.cor)[1]), y = 0:(dim(all.cor)[2]), z = 1-(all.cor^2), main = paste("Chr", ch), zlim = c(0, 1), xlim = c(0,dim(all.cor)[1]), ylim = c(0, dim(all.cor)[2]), xlab = "", ylab = "")
				# quartz();image(all.recomb)
				chr.block.locale <- which(sapply(strsplit(names(link.blocks), "_"), function(x) x[1]) == paste("Chr", ch, sep = ""))
				for(i in 1:length(chr.block.locale)){
					block.diag <- which(chr.markers %in% link.blocks[[chr.block.locale[i]]])-1
					block.diag <- c(block.diag, (max(block.diag)+1))
					lines(x = block.diag, y = block.diag, lwd = 2)
					lines(x = block.diag, y = rep(max(block.diag), length(block.diag)), lwd = 3)
					lines(x = rep(min(block.diag), length(block.diag)), y = block.diag, lwd = 3)
					}
				}
			}
		}		
	if(plot.blocks){dev.off()}
	
	if(collapse.linked.markers){
		data.obj$linkage.blocks.collapsed <- link.blocks
		}else{
		data.obj$linkage.blocks.full <- link.blocks
		}

	return(data.obj)
	
}
