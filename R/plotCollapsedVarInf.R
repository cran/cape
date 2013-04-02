plotCollapsedVarInf <-
function(data.obj, expand.labels = FALSE, all.markers = FALSE){
	
	adj.mat <- data.obj$collapsed.net
	if(is.null(adj.mat)){
		stop("This function operates on the collapsed network. collapse.net() must be run first.")
		}
	
	
	blocks <- data.obj$linkage.blocks.collapsed
	

	#replace marker numbers with names
	get.marker.names <- function(block){
		marker.locale <- which(colnames(data.obj$geno) %in% as.character(block))
		marker.names <- data.obj$marker.names[marker.locale]
		return(marker.names)
		}


	if(!all.markers){	

		adj.mat[which(adj.mat == 0)] <- NA
		
		if(expand.labels){
			named.blocks <- lapply(blocks, get.marker.names)			
			marker.names <- sapply(named.blocks, function(x) paste(x, collapse = ", "))
			pheno.names <- names(data.obj$max.var.to.pheno.influence)
			rownames(adj.mat) <- marker.names
			colnames(adj.mat) <- c(marker.names, pheno.names)
			}

		myImagePlot(adj.mat, min.x = (max(abs(adj.mat), na.rm = TRUE)*-1), max.x = max(abs(adj.mat), na.rm = TRUE), main = "Condensed Variant Influences", xlab = "Target", ylab = "Source")

		}else{ #if we are including all markers, we need to expand the adjacency matrix
			all.markers <- colnames(data.obj$geno)
			pheno.names <- names(data.obj$max.var.to.pheno.influence)
			expanded.adj.mat <- matrix(0, ncol = length(all.markers), nrow = length(all.markers))
			colnames(expanded.adj.mat) <- rownames(expanded.adj.mat) <- all.markers
			expanded.pheno.mat <- matrix(0, nrow = length(all.markers), ncol = length(pheno.names))
			colnames(expanded.pheno.mat) <- pheno.names; rownames(expanded.pheno.mat) <- all.markers
			expanded.adj.mat <- cbind(expanded.adj.mat, expanded.pheno.mat)
					
			#go through the blocks. Replace markers in blocks with the block number
			extra.markers.in.blocks <- sapply(blocks, length)
			big.block.locale <- which(extra.markers.in.blocks > 1)
			if(length(big.block.locale) > 0){
				for(i in 1:length(big.block.locale)){
					block.markers <- blocks[[big.block.locale[i]]]
					marker.locale <- which(colnames(expanded.adj.mat) %in% block.markers)
					col.to.remove <- marker.locale[-1]
					expanded.adj.mat <- expanded.adj.mat[,-col.to.remove]
					expanded.adj.mat <- expanded.adj.mat[-col.to.remove,]
					colnames(expanded.adj.mat)[marker.locale[1]] <- names(big.block.locale)[i]
					rownames(expanded.adj.mat)[marker.locale[1]] <- names(big.block.locale)[i]
					}
				}
			
							
			small.block.locale <- which(extra.markers.in.blocks == 1)		
			small.block.markers <- sapply(blocks, function(x) x[[1]])[small.block.locale]
			exp.block.locale <- which(rownames(expanded.adj.mat) %in% small.block.markers)
			rownames(expanded.adj.mat)[exp.block.locale] <- colnames(expanded.adj.mat)[exp.block.locale] <- names(small.block.markers)
			
			get.marker.name <- function(marker.number){
				marker.locale <- which(colnames(data.obj$geno) == marker.number)
				return(data.obj$marker.names[marker.locale])
				}
				
			#rename the rest of the markers with their names
			marker.rows <- rownames(expanded.adj.mat)
			block.locale <- grep("Block", marker.rows)
			not.block <- setdiff(marker.rows, marker.rows[block.locale])
			not.block.locale <- which(marker.rows %in% not.block)
			row.marker.names <- apply(matrix(not.block, nrow = 1), 2, get.marker.name)
			marker.rows[not.block.locale] <- row.marker.names

			rownames(expanded.adj.mat) <- marker.rows
			colnames(expanded.adj.mat)[1:length(marker.rows)] <- marker.rows

			
			sig.locale <- which(adj.mat != 0, arr.ind = TRUE)
			if(length(sig.locale) > 0){
				for(i in 1:length(sig.locale[,1])){
					location <- sig.locale[i,]
					adj.rowname <- rownames(adj.mat)[location[1]]
					adj.colname <- colnames(adj.mat)[location[2]]
					expanded.adj.mat[adj.rowname, adj.colname] <- adj.mat[location[1], location[2]]
					}
				}
						
		if(expand.labels){
			block.locale <- which(colnames(expanded.adj.mat) %in% names(blocks))
			if(length(block.locale) > 0){
				for(i in 1:length(block.locale)){
					block.name <- colnames(expanded.adj.mat)[block.locale[i]]
					block.markers <- blocks[[block.name]]
					block.marker.locale <- which(colnames(data.obj$geno) %in% block.markers)
					block.marker.name <- paste(data.obj$marker.names[block.marker.locale], collapse = ", ")
					colnames(expanded.adj.mat)[block.locale[i]] <- block.marker.name
					rownames(expanded.adj.mat)[block.locale[i]] <- block.marker.name
					}
				}	
			}
			
		expanded.adj.mat[which(expanded.adj.mat == 0)] <- NA
		myImagePlot(expanded.adj.mat, min.x = (max(abs(expanded.adj.mat), na.rm = TRUE)*-1), max.x = max(abs(expanded.adj.mat), na.rm = TRUE), main = "Condensed Variant Influences", xlab = "Target", ylab = "Source")

			
		}
	
	
	}
