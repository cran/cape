get.covar <-
function(data.obj, covar.thresh = NULL){
	
	oneD <- data.obj$singlescan.results

	if(is.null(covar.thresh)){
		covar.thresh <- data.obj$covar.thresh
		}else{
			data.obj$covar.thresh <- covar.thresh
			}
	
	covar.flags <- matrix(0, nrow = dim(oneD[[1]])[1], ncol = length(oneD))
	rownames(covar.flags) <- rownames(oneD[[1]])
	colnames(covar.flags) <- names(oneD)
		
	for(i in 1:dim(covar.flags)[2]){
		covar.flags[which(oneD[[i]][,"t.stat"] >= covar.thresh),i] <- 1
		}	

	data.obj$covar.flags <- covar.flags
	return(data.obj)
}
