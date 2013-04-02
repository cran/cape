set.covar <-
function(data.obj, pheno, markers, is.covar = TRUE, plot.covar = TRUE){
	
	
	covar.flags <- data.obj$covar.flags
	
	if(is.character(markers[1])){
		row.locale <- which(rownames(covar.flags) %in% markers)
		}else{
		row.locale <- markers
		}

	if(length(row.locale) < length(markers)){
		if(is.character(markers)[1]){
			didnt.find <- setdiff(markers, rownames(covar.flags)[row.locale])
			}else{
			didnt.find <- setdiff(markers, c(1:dim(covar.flags)[1])[row.locale])		
			}
		cat("\nI couldn't find the following markers:\n")
		cat(didnt.find, sep = "\n")
		return(data.obj)
		}


	if(is.character(pheno)[1]){
		col.locale <- which(colnames(covar.flags) %in% pheno)		
		}else{
		col.locale <- pheno	
		}

	if(length(col.locale) < length(pheno)){
		if(is.character(pheno)[1]){
			didnt.find <- setdiff(pheno, colnames(covar.flags)[col.locale])
			}else{
			didnt.find <- setdiff(pheno, c(1:dim(covar.flags)[2])[col.locale])	
			}
		cat("\nI couldn't find the following phenotypes:\n")
		cat(didnt.find, sep = "\n")
		return(data.obj)	
		}


	if(is.covar){
		covar.flags[row.locale, col.locale] <- 1
		}else{
		covar.flags[row.locale, col.locale] <- 0	
		}
		
	

	data.obj$covar.flags <- covar.flags

	if(plot.covar){
		plotSinglescan(data.obj, mark.covar = TRUE)
		}


	return(data.obj)

}
