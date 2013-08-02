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
	if(is.null(show.labels)){show.labels <- TRUE}
	additional.arguments$show.labels <- NULL
	chromosome.coordinates <- additional.arguments$chromosome.coordinates
	additional.arguments$chromosome.coordinates <- NULL
	chr.names <- additional.arguments$chr.names
	additional.arguments$chr.names <- NULL
	
	mark.coords <- additional.arguments$mark.coords
	mark.col <- additional.arguments$mark.col
	additional.arguments$mark.coords <- NULL
	additional.arguments$mark.col <- NULL


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
	ColorRamp <- rgb( seq(0,1,length=256),  # Red
	                   seq(0,1,length=256),  # Green
	                   seq(1,0,length=256))  # Blue

	ColorLevels <- seq(min.x, max.x, length=length(ColorRamp))

	
	#=====================================
	#plot the y axis label
	par(mar = c(0,0,5,0))
	plot.new()
	plot.window(xlim = c(0,1), ylim = c(0,1))
	if(show.labels){
		text(x = 0.3, y = 0.5, ylab, srt = 90)
		}else{
		text(x = 0.85, y = 0.5, ylab, srt = 90)	
		}
	
	if(show.labels){
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

	chr.cols <- rep(c("darkgray", "white"), ceiling(length(chromosome.coordinates)/2))
	y.chr.coord <- length(yLabels) - chromosome.coordinates + 1
	poly.perc = 0.03; label.size <- 0.9
	plot.dim <- par("usr")
	poly.width <- plot.dim[2]*poly.perc
	poly.max <- plot.dim[1]; poly.min <- poly.max-poly.width; poly.mid <- mean(c(poly.min, poly.max))


	if(!is.null(chromosome.coordinates)){
		par(xpd = TRUE)
		for(i in 1:(length(chromosome.coordinates)-1)){
				if(dim(x)[2] >= max(chromosome.coordinates)+1){
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
	
	# Color Scale
	par(mar = c(3,2.5,5,2))
	image(1, ColorLevels,
	matrix(data=ColorLevels, ncol=length(ColorLevels),nrow=1), col=ColorRamp, xlab="",ylab="",xaxt="n", cex.axis = 2)
	
	#plot the x axis labels	
	par(mar = c(0,0,0,2))
	plot.new()
	plot.window(xlim = c(0,1), ylim = c(0,1))
	if(show.labels){
		text(x = 0.5, y = 0.5, xlab)
		}else{
		text(x = 0.5, y = 0.85, xlab)
		}
}
