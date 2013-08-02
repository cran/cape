set.pairscan.thresh <-
function(data.obj, pairscan.thresh){
	
	
	if(!is.numeric(pairscan.thresh)){
		stop("pairscan.thresh must be a numeric value.")
		}
	
	data.obj$pairscan.thresh <- pairscan.thresh
	
	#change the value in alpha for pairs so it doesn't
	#continue to say the original alpha value
	#in the singlescan plots
	
	data.obj$alpha.for.pairs <- paste("t =", pairscan.thresh)
	
	return(data.obj)	
}
