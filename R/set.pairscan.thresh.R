set.pairscan.thresh <-
function(data.obj, pairscan.thresh){
	
	if(!is.numeric(pairscan.thresh)){
		stop("pairscan.thresh must be a numeric value.")
		}
	
	data.obj$pairscan.thresh <- pairscan.thresh
	return(data.obj)	
}
