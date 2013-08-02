### R code from vignette source 'cape.Rnw'

###################################################
### code chunk number 1: cape.Rnw:98-100
###################################################
library(cape)
data(obesity.cross)


###################################################
### code chunk number 2: cape.Rnw:131-132
###################################################
str(obesity.cross)


###################################################
### code chunk number 3: cape.Rnw:168-170
###################################################
obesity.cross <- select.pheno(obesity.cross, 
phenotypes = c("body_weight", "glucose", "insulin", "mom"))


###################################################
### code chunk number 4: cape.Rnw:179-180
###################################################
obesity.cross <- create.covar(obesity.cross, "mom")


###################################################
### code chunk number 5: cape.Rnw:205-209
###################################################
layout(matrix(c(1:3), nrow = 1))
hist.data <- apply(matrix(c(1:dim(obesity.cross$pheno)[2]), nrow = 1), 2, function(x) 
hist(obesity.cross$pheno[,x], main = colnames(obesity.cross$pheno)[x], 
xlab = colnames(obesity.cross$pheno)[x]))


###################################################
### code chunk number 6: cape.Rnw:221-233
###################################################
layout(matrix(c(1:3), nrow = 1))
i <- 1
while(i < dim(obesity.cross$pheno)[2]){
	j <- i + 1
	while(j <= dim(obesity.cross$pheno)[2]){
		qqplot(obesity.cross$pheno[,i], obesity.cross$pheno[,j], 
		xlab = colnames(obesity.cross$pheno)[i], ylab = 
		colnames(obesity.cross$pheno)[j], cex.lab = 1.5)
		j <- j + 1
	}
	i = i + 1
}


###################################################
### code chunk number 7: cape.Rnw:246-247
###################################################
obesity.cross <- norm.pheno(obesity.cross, mean.center = TRUE)


###################################################
### code chunk number 8: cape.Rnw:254-258
###################################################
layout(matrix(c(1:3), nrow = 1))
hist.data <- apply(matrix(c(1:dim(obesity.cross$pheno)[2]), nrow = 1), 2, function(x) 
hist(obesity.cross$pheno[,x], main = colnames(obesity.cross$pheno)[x], 
xlab = colnames(obesity.cross$pheno)[x]))


###################################################
### code chunk number 9: cape.Rnw:268-280
###################################################
layout(matrix(c(1:3), nrow = 1))
i <- 1
while(i < dim(obesity.cross$pheno)[2]){
	j <- i + 1
	while(j <= dim(obesity.cross$pheno)[2]){
		qqplot(obesity.cross$pheno[,i], obesity.cross$pheno[,j], 
		xlab = colnames(obesity.cross$pheno)[i], ylab = 
		colnames(obesity.cross$pheno)[j], cex.lab = 1.5)
		j <- j + 1
	}
	i = i + 1
}


###################################################
### code chunk number 10: cape.Rnw:313-327
###################################################
layout(matrix(c(1:3), nrow = 1))
i <- 1
while(i < dim(obesity.cross$pheno)[2]){
	j <- i + 1
	while(j <= dim(obesity.cross$pheno)[2]){
		plot(obesity.cross$pheno[,i], obesity.cross$pheno[,j], 
		xlab = colnames(obesity.cross$pheno)[i], ylab = 
		colnames(obesity.cross$pheno)[j], main = paste("r =", 
		round(cor(obesity.cross$pheno[,i], obesity.cross$pheno[,j], 
		use = "complete.obs"), 2)), cex.lab = 1.2)
		j <- j + 1
	}
	i = i + 1
}


###################################################
### code chunk number 11: cape.Rnw:362-364
###################################################
obesity.cross <- get.eigentraits(obesity.cross, scale.pheno = FALSE,
normalize.pheno = FALSE)


###################################################
### code chunk number 12: cape.Rnw:372-373
###################################################
plotSVD(obesity.cross, orientation = "vertical")


###################################################
### code chunk number 13: cape.Rnw:413-414
###################################################
obesity.cross <- select.eigentraits(obesity.cross, traits.which = c(1,2))


###################################################
### code chunk number 14: cape.Rnw:448-450
###################################################
obesity.cross <- singlescan(obesity.cross, n.perm = 100, covar = "mom", 
scan.what = "eigentraits", alpha = c(0.01, 0.05), verbose = FALSE)


###################################################
### code chunk number 15: cape.Rnw:484-485
###################################################
plotSinglescan(obesity.cross, mark.chr = TRUE, mark.covar = FALSE)


###################################################
### code chunk number 16: cape.Rnw:552-554
###################################################
obesity.cross <- select.markers.for.pairscan(obesity.cross, 
use.pairs.threshold = FALSE)


###################################################
### code chunk number 17: cape.Rnw:639-641
###################################################
plotSinglescan(obesity.cross, mark.chr = TRUE, show.rejected.markers = TRUE, 
standardized = TRUE)


