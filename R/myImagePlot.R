myImagePlot <-
function(x,...){
	# print(dim(x))

	#build the argument list from additional arguments added to
	#the function
	additional.arguments <- list(...)

	#=================================================================
	#There are some special additional arguments that this
	#function can take in
	#if there are xlab and ylab arguments here
	#pull them out first before we override them
	xlab <- additional.arguments$xlab
	ylab <- additional.arguments$ylab
	show.labels <- additional.arguments$show.labels
	additional.arguments$show.labels <- NULL
	if(is.null(show.labels)){show.labels <- TRUE}

	show.pheno.labels <- additional.arguments$show.pheno.labels
	additional.arguments$show.pheno.labels <- NULL
	if(is.null(show.pheno.labels)){show.pheno.labels <- TRUE}
		
	chromosome.coordinates <- additional.arguments$chromosome.coordinates
	additional.arguments$chromosome.coordinates <- NULL
	chr.names <- additional.arguments$chr.names
	additional.arguments$chr.names <- NULL
	
	mark.coords <- additional.arguments$mark.coords
	mark.col <- additional.arguments$mark.col
	additional.arguments$mark.coords <- NULL
	additional.arguments$mark.col <- NULL

	extra.col.mat <- additional.arguments$extra.col.mat
	additional.arguments$extra.col.mat <- NULL

	#if there are min.x and max.x argments
	#pull these out too.
	min.x <- additional.arguments$min.x
	max.x <- additional.arguments$max.x
	#=================================================================
	
	
	#if they aren't specified, use the
	#x matrix to specify them
	if(is.null(min.x)){
		min.x <- min(x, na.rm = TRUE)
		}else{ #otherwise remove it from the argument list, so it doesn't throw a warning when we use it in image
			additional.arguments$min.x <- NULL
			}
	if(is.null(max.x)){
		max.x <- max(x, na.rm = TRUE)
		}else{
			additional.arguments$max.x <- NULL
			}


	yLabels <- rownames(x)
	xLabels <- colnames(x)

	layout.mat <- matrix(c(1:3, 4, 4, 4), nrow=2, ncol=3, byrow = TRUE)
	layout(layout.mat, widths=c(0.75,4,0.75), heights = c(1,0.1))
	# layout.show(4);return()
	
	# Red and green range from 0 to 1 while Blue ranges from 1 to 0
	# ColorRamp <- rgb(seq(0,1,length=256),  # Red
	                   # seq(0,1,length=256),  # Green
	                   # seq(1,0,length=256))  # Blue
	
	# ColorLevels <- seq(min.x, max.x, length=length(ColorRamp))
	ColorLevels <- seq(min.x, max.x, length=256)

	#center the palette on 0
	# mypal.neg <- colorRampPalette(c("dodgerblue3", "aliceblue"))
	# mypal.pos <- colorRampPalette(c("lavenderblush", "red"))

	#spectral colors
	# mypal.pos <- colorRampPalette(rev(c("#9E0142", "#D53E4F", "#F46D43", "#FDAE61", "#FEE08B", "#FFFFBF")))
	# mypal.neg <- colorRampPalette(rev(c("#FFFFBF", "#E6F598", "#ABDDA4", "#66C2A5", "#3288BD", "#5E4FA2")))
	
	#purple-green
	# mypal.pos <- colorRampPalette(rev(c("#40004B", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#F7F7F7")))
	# mypal.neg <- colorRampPalette(rev(c("#F7F7F7", "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B")))

	#purple-green2
	mypal.pos <- colorRampPalette(c("#f7f7f7", "#af8dc3"))
	mypal.neg <- colorRampPalette(c("#7fbf7b", "#f7f7f7"))


	#pink-green
	# mypal.pos <- colorRampPalette(rev(c("#8E0152", "#C51B7D", "#DE77AE", "#F1B6DA", "#FDE0EF", "#F7F7F7")))
	# mypal.neg <- colorRampPalette(rev(c("#F7F7F7", "#E6F5D0", "#B8E186", "#7FBC41", "#4D9221", "#276419")))
	# mypal.pos <- colorRampPalette(c("#f7f7f7","#e9a3c9"))
	# mypal.neg <- colorRampPalette(c("#a1d76a","#f7f7f7"))

	
	#purple-orange
	# mypal.pos <- colorRampPalette(rev(c("#7F3B08", "#B35806", "#E08214", "#FDB863", "#FEE0B6", "#F7F7F7")))
	# mypal.neg <- colorRampPalette(rev(c("#F7F7F7", "#D8DAEB", "#B2ABD2", "#8073AC", "#542788", "#2D004B")))  

	#purple-orange2
	# mypal.pos <- colorRampPalette(c("#f7f7f7", "#f1a340"))
	# mypal.neg <- colorRampPalette(c("#998ec3", "#f7f7f7"))  

	#red-blue
 	# mypal.pos <- colorRampPalette(c("#f7f7f7", "#ef8a62"))
 	# mypal.neg <- colorRampPalette(c("#67a9cf" ,"#f7f7f7"))  

	#red-yellow-blue
 	# mypal.pos <- colorRampPalette(rev(c("#A50026", "#D73027", "#F46D43", "#FDAE61", "#FEE090", "#FFFFBF")))
 	# mypal.neg <- colorRampPalette(rev(c("#FFFFBF", "#E0F3F8", "#ABD9E9", "#74ADD1", "#4575B4", "#313695")))  
	
	#green-gray-yellow
	# mypal.pos <- colorRampPalette(c("midnightblue", "springgreen"))
	# mypal.neg <- colorRampPalette(c("yellow", "midnightblue"))

	# mypal.pos <- colorRampPalette(c("gray", "red"))
	# mypal.neg <- colorRampPalette(c("blue", "gray"))

	#blue brown
	# mypal.pos <- colorRampPalette(c("#FFFFFF", "#3794bf"))
	# mypal.neg <- colorRampPalette(c("#df8640", "#FFFFFF"))
	# mypal.pos <- colorRampPalette(rev(c("#a6611a", "#dfc27d", "#f5f5f5")))
	# mypal.neg <- colorRampPalette(rev(c("#f5f5f5", "#80cdc1", "#018571")))
	

	#jet colors
	# mypal.neg <- colorRampPalette(c("blue", "#007FFF", "cyan","#7FFF7F"))
	# mypal.pos <- colorRampPalette(c("#7FFF7F", "yellow", "#FF7F00", "red"))

	# mypal.pos <- colorRampPalette(c("white", "#fc8d59"))
	# mypal.neg <- colorRampPalette(c("#91bfdb", "white"))

	
	ColorRamp <- c(mypal.neg(length(which(ColorLevels < 0))), mypal.pos(length(which(ColorLevels >= 0))))
	
	#=====================================
	#plot the y axis label
	par(mar = c(0,0,5,0))
	plot.new()
	plot.window(xlim = c(0,1), ylim = c(0,1))
	if(show.labels){
		text(x = 0.3, y = 0.5, ylab, srt = 90, cex = 2)
		}else{
		text(x = 0.85, y = 0.5, ylab, srt = 90, cex = 2)	
		}
	
	if(show.labels || show.pheno.labels){
		par(mar = c(5,3,5,2))
		}else{
		par(mar = c(3,3,5,2))
		}

	#add the default arguments to the argument list
	additional.arguments$x <- 1:length(xLabels)
	additional.arguments$y <- 1:length(yLabels)
	additional.arguments$z <- rotate.mat(x)
	# additional.arguments$z <- x
	additional.arguments$col = ColorRamp
	additional.arguments$xlab <- ""
	additional.arguments$ylab = ""
	additional.arguments$axes = FALSE
	additional.arguments$zlim <- c(min.x, max.x)
	# additional.arguments$useRaster <- TRUE
	do.call(image, additional.arguments)

	#add the extra colors if we are highlighting particular cells
	if(!is.null(extra.col.mat)){
		u_col <- unique(as.vector(extra.col.mat[which(!is.na(extra.col.mat))]))
		if(length(u_col) > 0){
			for(i in 1:length(u_col)){
				col.locale <- which(rotate.mat(extra.col.mat) == u_col[i], arr.ind = TRUE)
				points(col.locale[,1], col.locale[,2], col = u_col[i], pch = 15, cex = 0.6)
				}
			}
		}


	chr.cols <- rep(c("darkgray", "white"), ceiling(length(chromosome.coordinates)/2))
	y.chr.coord <- length(yLabels) - chromosome.coordinates + 1
	poly.perc = 0.03; label.size <- 0.9
	plot.dim <- par("usr")
	poly.width <- plot.dim[2]*poly.perc
	poly.max <- plot.dim[1]; poly.min <- poly.max-poly.width; poly.mid <- mean(c(poly.min, poly.max))


	if(!is.null(chromosome.coordinates)){
		par(xpd = TRUE)
		for(i in 1:(length(chromosome.coordinates)-1)){
				if(dim(x)[2]+1 >= max(chromosome.coordinates)){
				polygon(x = c(chromosome.coordinates[i], chromosome.coordinates[i+1], chromosome.coordinates[i+1], chromosome.coordinates[i]), y = c(poly.min, poly.min, poly.max, poly.max), col = chr.cols[i])
				text(x = mean(c(chromosome.coordinates[i], chromosome.coordinates[i+1])), y = poly.mid, cex = label.size, labels = chr.names[i])
				}

				polygon(y = c(y.chr.coord[i], y.chr.coord[i+1], y.chr.coord[i+1], y.chr.coord[i]), x = c(poly.min, poly.min, poly.max, poly.max), col = chr.cols[i])
				text(y = mean(c(y.chr.coord[i], y.chr.coord[i+1])), x = poly.mid, cex = label.size, labels = chr.names[i], srt = 90)
			
			}
			
		par(xpd = FALSE)
		}

	if(!is.null(mark.coords)){
		points(mark.coords[,1], mark.coords[,2], col = mark.col, pch = 16)
		}
	
	
	if(show.labels){
		par(xpd = TRUE)
		if(is.null(chromosome.coordinates)){
			axis(LEFT<-2, at=1:length(yLabels), labels=rev(yLabels), las= HORIZONTAL<-1,cex.axis=0.7, adj = 1, tick = TRUE)
			axis(BELOW<-1, at=1:length(xLabels), labels=FALSE, cex.axis=0.7, tick = TRUE)
			text(x = 1:length(xLabels), y = poly.min*1.2, labels = xLabels, srt = 90, cex = 0.7, adj = 1)
			}else{
			axis(LEFT<-2, at=1:length(yLabels), labels=rev(yLabels), las= HORIZONTAL<-1,cex.axis=0.7, adj = 1, tick = FALSE, line = 0.8)
			axis(BELOW<-1, at=1:length(xLabels), labels=FALSE, cex.axis=0.7, tick = FALSE)
			text(x = 1:length(xLabels), y = poly.min*1.2, labels = xLabels, srt = 90, cex = 0.7, adj = 1)
			}
		par(xpd = FALSE)
		}else{
		axis(LEFT<-2, at=1:length(yLabels), labels=FALSE, las= HORIZONTAL<-1,cex.axis=0.7, adj = 1, tick = TRUE, lwd.ticks = 0)
		axis(BELOW<-1, at=1:length(xLabels), labels=FALSE, cex.axis=0.7, tick = TRUE, lwd.ticks = 0)
		}
	
		if(show.pheno.labels & !show.labels){
			only.pheno.labels <- xLabels
			marker.locale <- which(yLabels %in% xLabels)
			only.pheno.labels[marker.locale] <- ""
			par(xpd = TRUE)
			text(x = 1:length(only.pheno.labels), y = mean(c(poly.min, poly.max)), labels = only.pheno.labels, srt = 90, cex = 1.5, adj = 1)
			par(xpd = FALSE)
			}

	#plot the x axis label close to the axis if we are not printing the labels
	if(!show.labels){
		par(xpd = TRUE)
		plot.height = plot.dim[2]-plot.dim[1]
		marker.locale <- which(yLabels %in% xLabels)
		text(x = median(marker.locale), y = poly.min-(plot.height*0.05), xlab, cex = 2)
		par(xpd = FALSE)
		}



	# Color Scale
	par(mar = c(3,2.5,5,2))
	image(1, ColorLevels, matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp, xlab="",ylab="",xaxt="n", cex.axis = 2)
		
	#plot the x axis labels in a different window if we are printing the marker labels

	par(mar = c(0,0,0,2))
	plot.new()
	plot.window(xlim = c(0,1), ylim = c(0,1))
	if(show.labels){
		par(xpd = TRUE)
		text(x = 0.5, y = 0.5, xlab, cex = 2)
		par(xpd = FALSE)
		}

}
