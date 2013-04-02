plotNetwork <-
function(data.obj, collapsed.net = TRUE, trait = NULL, chr = NULL){

	require(igraph)	
	
	iArrows <- igraph:::igraph.Arrows
			
	all.chr <- data.obj$chromosome
	all.pos <- data.obj$marker.location

	
	if(collapsed.net){
		adj.mat <- data.obj$collapsed.net
		}else{
		adj.mat <- data.obj$full.net	
		}
		
		
		if(is.null(adj.mat)){
			stop("get.network() must be run before plotting the collapsed network.")
			}
		
		if(collapsed.net){
			blocks <- data.obj$linkage.blocks.collapsed
			}else{
			blocks <- data.obj$linkage.blocks.full	
			}
			
		all.markers <- as.vector(unlist(blocks))
		
		
		#assign a chromosome and relative position to each block
		get.chr.pos <- function(block){
			marker.locale <- which(colnames(data.obj$geno) %in% block)
			chr <- unique(all.chr[marker.locale])
			pos <- mean(all.pos[marker.locale])
			max.pos <- max(all.pos[all.chr == chr])
			return(c(chr, pos/max.pos))
			}
				
		chr.pos <- t(sapply(blocks, get.chr.pos))
		colnames(chr.pos) <- c("chromosome", "position")

		if(is.null(chr)){
			chr <- unique(all.chr)
			}
			
		#calculate beginning and end x coordinates for each chromosome
		num.chr <- length(chr)
		chr.x <- matrix(NA, ncol = 2, nrow = num.chr)
		rownames(chr.x) <- chr
		chr.x[,1] <- (0:(num.chr-1))+(num.chr*0.005)
		chr.x[,2] <- (1:num.chr)-(num.chr*0.005)
		
		
		#and each phenotype
		if(is.null(trait)){
			pheno <- names(data.obj$max.var.to.pheno.influence)
			}else{
			pheno <- trait
			trait.locale <- which(trait %in% names(data.obj$max.var.to.pheno.influence))
			if(length(trait.locale) < length(trait)){
				not.found <- which(!trait %in% names(data.obj$max.var.to.pheno.influence))
				message("I couldn't find the following traits:")
				cat(trait[not.found], sep = "\n")
				return()
				}
			}
			
		
		#if we need to filter chr.pos and adj.mat to include
		#only the chromomsomes and phenotypes we are including
		chr.pos <- chr.pos[which(chr.pos[,"chromosome"] %in% chr),]
		all.block.pheno <- c(rownames(chr.pos), pheno)
		adj.mat <- adj.mat[,colnames(adj.mat) %in% all.block.pheno]
		adj.mat <- adj.mat[rownames(adj.mat) %in% rownames(chr.pos),]
		
		#get the absolute position of each marker
		get.abs.mrk <- function(chr.pos.row){
			chr.locale <- which(rownames(chr.x) == chr.pos.row[1])
			mrk.pos <- ((chr.x[chr.locale,2] - chr.x[chr.locale,1])*as.numeric(chr.pos.row[2])) + chr.x[chr.locale,1]
			return(as.numeric(mrk.pos))
			}
		
		ph.x <- apply(chr.pos, 1, get.abs.mrk)
		
		#if there are covariates, expand the covariate "chromosome"
		#to give each covariate its own segment
		covar.exists <- which(chr.pos[,1] == 0)
		covar.locale <- which(rownames(chr.x) == 0)
		if(length(covar.locale) > 0){
			covar.min <- chr.x[covar.locale,1]
			covar.max <- chr.x[covar.locale,2]
			new.cov.x <- matrix(NA, ncol = 2, nrow = length(which(all.chr == 0)))
			region.centers <- segment.region(covar.min, covar.max, length(which(all.chr == 0)), "center") 
			new.cov.x[,1] <- region.centers - ((covar.max-covar.min)/100)
			new.cov.x[,2] <- region.centers + ((covar.max-covar.min)/100)
			rownames(new.cov.x) <- rep(0, dim(new.cov.x)[1])
			chr.x <- chr.x[which(rownames(chr.x) != 0),]
			chr.x <- rbind(chr.x, new.cov.x)
			
			#also adjust the positions of the covariates in ph.x
			covar.blocks <- rownames(chr.pos)[which(chr.pos[,1] == 0)]
			if(length(covar.blocks) > 0){
				sig.covar <- NULL
				for(i in 1:length(covar.blocks)){
					sig.covar <- blocks[[covar.blocks[i]]]
					}
				covar.pos <- which(colnames(data.obj$geno)[which(all.chr == 0)] %in% sig.covar)
				cov.centers <- segment.region(covar.min, covar.max, length(which(all.chr == 0)), "center")
				cov.locale <- which(chr.pos[,1] == 0)
				ph.x[cov.locale] <- cov.centers[covar.pos]
				}
			}
		
		
		#take out the covariate label and add it in at the end
		chr.labels <- chr
		chr.labels[which(chr.labels == 0)] <- ""

		dev.new(width = 12, height = 7)
		plot.new()
		plot.window(xlim = c(0,length(chr)), ylim = c(0,1))
		for(i in 1:length(chr.x[,1])){
			segments(chr.x[i,1], 0.5, chr.x[i,2], 0.5, lwd = 5)
			text(x = mean(chr.x[i,]), y = 1, chr.labels[i], cex = 1.5)
			}
		if(length(covar.locale) > 0){
			text(x = mean(new.cov.x), y = 1, "Cov.", cex = 1.5)
			}
			
		par(xpd = TRUE)
		num.pheno <- length(pheno)
		ph.y <- 0.45
		for(ph in 1:length(pheno)){
			edge.col <- rep("gray", length(ph.x))
			sig.locale <- which(adj.mat[,pheno[ph]] != 0)
			if(length(sig.locale) > 0){
				edge.col[which(adj.mat[,pheno[ph]] > 0)] <- rgb(0, 201/256, 87/256)
				edge.col[which(adj.mat[,pheno[ph]] < 0)] <- "red"
				}
			text(x = -0.8, y = ph.y, pheno[ph], cex = 0.8)
			points(x = ph.x, y = rep(ph.y, length(ph.x)), col = edge.col, pch = "|", cex = 1.5)
			ph.y <- ph.y - 0.05
			}
		

		#add arrows for the variant interactions
		just.m <- adj.mat[,-which(colnames(adj.mat) %in% pheno)]
		for(i in 1:dim(just.m)[1]){
			sig.locale <- which(just.m[i,] != 0)
			
			if(length(sig.locale) > 0){
				edge.col <- rep(NA, length(sig.locale))
				edge.col[which(just.m[i,sig.locale] > 0)] <- rgb(0, 201/256, 87/256)
				edge.col[which(just.m[i,sig.locale] < 0)] <- "red"
				
				target.spec.pos <- ph.x[colnames(just.m)[sig.locale]]
				source.spec.pos <- rep(ph.x[rownames(just.m)[i]], length(sig.locale))
				
				curve.dir <- sign(target.spec.pos-source.spec.pos)
				curve.mag <- 1.3/length(chr)
				iArrows(x1 = source.spec.pos, y1 = rep(0.5, length(sig.locale)), x2 = target.spec.pos, y2 = rep(0.5, length(sig.locale)), h.lwd=1, sh.lwd=2, sh.col=edge.col, curve = curve.mag*curve.dir, width=1, size=0.7)
				}
			}
			
			
			if(collapsed.net){
				mtext("Linkage Block Influences", cex = 1.2, line = 2)
				}else{
				mtext("Variant Influences", cex = 1.2, line = 2)	
				}	

	
	}
