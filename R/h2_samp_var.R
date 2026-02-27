#' Calculates expected sampling variance of a given $h^2$ from the A matrix of a pedigree, based on the method of Visscher and Goddard 2015 Genetics

#'
#' @param Ped A pedigree
#' @param h2 vector of heritabilities for which the sampling variance is calculated.
#' @param plot logical. If TRUE, a plot of the sampling variance as a function of h2 is produced (only if length(h2)>1).
#'
#' @return
#'   Returns a vector of sampling variances for each given h2.
#'
#' @examples
#'\donttest{
#' data(gryphons)
#' pedigree <- gryphons[, 1:3]
#'
#' h2_samp_var(pedigree,h2=0.3)
#' }
#' @export
#'

h2_samp_var <- function(Ped, h2 = 0.3, plot = FALSE) {
	Ped <- fix_ped(Ped)
	A <- nadiv::makeA(Ped)
	lambda <- eigen(A, symmetric = TRUE, only.values = TRUE)$values
	N <- nrow(Ped)
	h2_samp_var <- sapply(h2, function(x) {
		a <- sum((lambda - 1)^2 / (1 + x * (lambda - 1))^2)
		b <- sum((lambda - 1) / (1 + x * (lambda - 1)))
		2 / (a - b^2 / N)
	})

	if (plot & length(h2) > 1) {
		plot(
			h2_samp_var ~ h2,
			xlab = "Heritability",
			ylab = "Expected Sampling Variance",
			main = "Expected Sampling Variance of Heritability based on Visscher & Goddard 2016",
			type = "b"
		)
	}
	names(h2_samp_var) <- h2
	return(h2_samp_var)
}
