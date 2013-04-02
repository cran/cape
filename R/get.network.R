get.network <-
function(data.obj, p.or.q = 0.05, collapse.linked.markers = TRUE, r2.thresh = 0.8, standardize = FALSE){
	
	if(!collapse.linked.markers){
		r2.thresh = 2
		}
	
	#get the linkage blocks based on the significant markers
	data.obj <- linkage.blocks(data.obj, p.or.q = p.or.q, collapse.linked.markers = collapse.linked.markers, r2.thresh = r2.thresh)
	
	if(collapse.linked.markers){
		blocks <- data.obj$linkage.blocks.collapsed
		}else{
		blocks <- data.obj$linkage.blocks.full	
		}
	
	if(length(blocks) == 1){
		stop("There is only one linkage block at this r2 threshold.")
		}
	
	#build a new network based on the block structure
	all.net.data <- data.obj$var.to.var.p.val
	
	var.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(all.net.data))))
	net.data <- all.net.data[which(as.numeric(all.net.data[,var.sig.col]) <= p.or.q),,drop = FALSE]
	pheno.tables <- data.obj$max.var.to.pheno.influence
	phenotypes <- names(pheno.tables)	

	adj.mat <- matrix(0, ncol = length(blocks), nrow = length(blocks))
	colnames(adj.mat) <- rownames(adj.mat) <- names(blocks)
	
	block.pairs <- pair.matrix(names(blocks))

	#for each pair of blocks
	get.adj.weight <- function(block.pair){
		#get all the markers in the two blocks
		all.markers1 <- blocks[[block.pair[1]]]
		all.markers2 <- blocks[[block.pair[2]]]
		
		#find the maximum weight between markers in the blocks
		#in both directions
		block1.source.locale <- which(net.data[,"Source"] %in% all.markers1)
		block2.target.locale <- which(net.data[,"Target"] %in% all.markers2)
		block1.to.block2 <- intersect(block1.source.locale, block2.target.locale)
		
		if(length(block1.to.block2) > 0){
			if(standardize){
				all.effects <- as.numeric(net.data[block1.to.block2,"Effect"])/as.numeric(net.data[block1.to.block2,"SE"])
				}else{
				all.effects <- as.numeric(net.data[block1.to.block2,"Effect"])	
				}
			adj.mat[block.pair[1], block.pair[2]] <- all.effects[which(abs(all.effects) == max(abs(all.effects)))][1]
			}
		

		block2.source.locale <- which(net.data[,"Source"] %in% all.markers2)
		block1.target.locale <- which(net.data[,"Target"] %in% all.markers1)
		block2.to.block1 <- intersect(block2.source.locale, block1.target.locale)
		
		if(length(block2.to.block1) > 0){
			if(standardize){
				all.effects <- as.numeric(net.data[block2.to.block1,"Effect"])/as.numeric(net.data[block2.to.block1,"SE"])
				}else{
				all.effects <- as.numeric(net.data[block2.to.block1,"Effect"])	
				}
			adj.mat[block.pair[2], block.pair[1]] <- all.effects[which(abs(all.effects) == max(abs(all.effects)))][1]
			}
			
		return(adj.mat)
		}
	
	for(i in 1:length(block.pairs[,1])){
		adj.mat <- get.adj.weight(block.pairs[i,])
		}
	
	
	#Now add the phenotypic effects continuing to use the maximum significant effect from each block
	pheno.mat <- matrix(0, nrow = length(blocks), ncol = length(phenotypes))
	colnames(pheno.mat) <- phenotypes
	rownames(pheno.mat) <- names(blocks)
	
	pheno.sig.col <- as.vector(na.omit(match(c("qval", "lfdr", "p.adjusted"), colnames(pheno.tables[[1]]))))
	get.block.inf <- function(block){
		all.markers <- blocks[[block]]
		for(i in 1:length(pheno.tables)){
			sig.inf <- pheno.tables[[i]][which(pheno.tables[[i]][,pheno.sig.col] <= p.or.q),,drop = FALSE]
			if(length(sig.inf) > 0){
				block.locale <- which(sig.inf[,"marker"] %in% all.markers)
				if(length(block.locale)){
					all.effects <- as.numeric(sig.inf[block.locale,"coef"])/as.numeric(sig.inf[block.locale,"se"])
					pheno.mat[block,i] <- all.effects[which(abs(all.effects) == max(abs(all.effects)))][1]
					}
				}
			}
		return(pheno.mat)	
		}

	for(i in 1:length(blocks)){
		pheno.mat <- get.block.inf(names(blocks)[i])
		}
	
	final.mat <- cbind(adj.mat, pheno.mat)

	if(collapse.linked.markers){
		data.obj$collapsed.net <- final.mat
		}else{
		data.obj$full.net <- final.mat	
		}
	
	return(data.obj)
	
	}
