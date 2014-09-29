
#' computes variograms for each of the coefficients
#'
#' @param geodat list of geodata objects. Each element of the list corresponds to a coefficient
#' @param estimator.type string or vector string specifying estimator type
#' @param max.dist scalar or vector specifying maximum distance to compute variogram
#' @param ... arguments passed to variog
variograms <- function( geodat, estimator.type="modulus", max.dist=0.6,...){
	n <- length(geodat)
	if(length(estimator.type)==1) estimator.type <- rep(estimator.type, n)
	if(length(max.dist)==1) max.dist <- rep(max.dist, n)
	
	res <- llply(1:n, 
		function(i, list, estimator, max.dist){
			variog(list[[i]], estimator.type=estimator[i], max.dist=max.dist[i], ...)
		}, list = geodat, estimator=estimator.type, max.dist=max.dist)
}