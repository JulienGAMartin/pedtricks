#' Calculates expected sampling variance of a given $h^2$ from the A matrix of a pedigree, based on the metric of Visscher and Goddard 2015 Genetics

#'
#' @param Ped A pedigree
#' @param h2 vector of heritabilities for which the sampling variance is calculated. 
#'
#' @return
#'   Returns a vector of sampling variances for each given h2.
#'
#' @examples
#'
#' data(gryphons)
#' pedigree <- gryphons[, 1:3]
#'
#' vg_samp_var(pedigree,h2=0.3)
#' 
#' @export
#'

vg_samp_var <- function(ped, h2){
	ped<-fix_ped(ped)
	A<-nadiv::makeA(ped)
	lambda <- eigen(A, symmetric=TRUE,only.values=TRUE)$values
	N <- nrow(ped)
	sapply(h2, function(x){
		a <- sum( (lambda-1)^2 / (1 + x*(lambda-1))^2 )
		b <- sum( (lambda-1) / (1 + x*(lambda-1)) )
		2/(a-b^2/N)
	})	
}
