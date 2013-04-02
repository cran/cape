linkage.blocks <-
function(data.obj, p.or.q = 0.05, collapse.linked.markers = TRUE, r2.thresh = 0.8){
	
	net.data <- data.obj$var.to.var.p.val
	pheno.net.data <- data.obj$max.var.to.pheno.influence
	data.obj$r2.thresh <- r2.thresh
	
	
	if(!collapse.linked.markers){
		r2.thresh <- 2
		}
	
	if(length(net.data) == 0){
		stop("calc.p() must be run to calculate variant-to-variant influences.")
		}
		
	
	#get the markers with significant influences on other markers
	#and on phenotypes
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(net.data))))
	sig.markers <- net.data[which(net.data[, var.sig.col] <= p.or.q), 1:2]
	
	if(length(sig.markers) == 0){
		stop("There are no significant blocks at this p.or.q value.")
		}

	sig.pheno <- vector(mode = "list", length = length(pheno.net.data))
	pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.net.data[[1]]))))
	for(ph in 1:length(pheno.net.data)){
		sig.pheno[[ph]] <- pheno.net.data[[ph]][which(pheno.net.data[[ph]][,pheno.sig.col] <= p.or.q),1]
		}

	if(length(c(sig.markers, unlist(sig.pheno))) == 0){
		stop("There are no significant markers at p = ", p.or.q)
		}


	#find all the markers with significant p values
	#and sort them
	u_markers <- unique(c(sig.markers, unlist(sig.pheno)))
	marker.locale <- match(u_markers, colnames(data.obj$geno))
	sorted.markers <- u_markers[order(marker.locale)]

	
	#take these out of the genotype matrix
	marker.geno <- data.obj$geno[,as.character(sorted.markers)]
	all.cor <- cor(marker.geno, use = "complete")^2
	#zero out the diagonal and the lower triangle
	all.cor[lower.tri(all.cor, diag = TRUE)] <- 0
	
	in.ld <- which(all.cor >= r2.thresh, arr.ind = TRUE)


	#start with a linkage block list that contains one marker per block
	link.block <- as.vector(sorted.markers, mode = "list")
	names(link.block) <- paste("Block", 1:length(link.block), sep = "")
	
	# name.blocks <- lapply(link.block, get.marker.names)
	
	#if there are no markers in LD, just return the list as is.
	if(length(in.ld) == 0){
		if(collapse.linked.markers){
			data.obj$linkage.blocks.collapsed <- link.block
			}else{
			data.obj$linkage.blocks.full <- link.block
			}
			return(data.obj)
		}
	
	#otherwise, go through the in.ld matrix and combine markers that are linked
	linked.markers <- list(in.ld[1,])
	
	if(dim(in.ld)[1] > 1){
		for(i in 2:length(in.ld[,1])){

			row.to.check <- in.ld[i,]
			#look for blocks that already contain markers we're looking at
			#The lapply function subtracts the length of the combined unique
			#vector from the length of the total vector. If the result is 
			#positive, there are common nodes between the new row and an
			#existing block
			common.nodes <- lapply(linked.markers, function(x) length(c(x,row.to.check))-length(unique(c(x,row.to.check))))
			shared.node.locale <- which(common.nodes > 0)

			#if there are no shared nodes, start a new block 
			if(length(shared.node.locale) == 0){
				linked.markers[[(length(linked.markers)+1)]] <- as.numeric(row.to.check)
				}

			#if we need to combine the new row with one other block
			if(length(shared.node.locale) == 1){
				linked.markers[[shared.node.locale]] <- as.numeric(unique(c(linked.markers[[shared.node.locale]], row.to.check)))
				}
		
			if(length(shared.node.locale) > 1){
				message("There is more than one common node")
				}

			}
		}

	#go through the linkage blocks and adjust the marker block list
	for(i in 1:length(linked.markers)){
		marker.names <- as.numeric(colnames(marker.geno)[sort(linked.markers[[i]])])
		# cat(i, marker.names, "\n")
		list.locale <- match(marker.names, link.block)
		#replace the first marker in the block with all markers
		link.block[[min(list.locale)]] <- marker.names
		#nullify the other blocks
		to.nullify <- list.locale[-1]
		for(j in 1:length(to.nullify)){
			# cat("\ttaking out", link.block[[to.nullify[j]]], "\n")
			link.block[[to.nullify[j]]] <- 0
			}
		}
	#remove all the blocks that have been converted to 0
	to.nullify <- which(link.block %in% 0)
	while(length(to.nullify) > 0){
		link.block[[to.nullify[1]]] <- NULL
		to.nullify <- which(link.block %in% 0)
		}
	
			
	# name.blocks <- lapply(link.block, get.marker.names)
	# names(name.blocks) <- paste("Block", 1:length(link.block), sep = "")
	
	
	if(collapse.linked.markers){
		data.obj$linkage.blocks.collapsed <- link.block
		}else{
		data.obj$linkage.blocks.full <- link.block
		}
	return(data.obj)
	
}
