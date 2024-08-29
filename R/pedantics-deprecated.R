#' Deprecated functions from package \pkg{pedantics}.
#'
#' The functions listed below are deprecated and are available for
#'   backward compatibility only. They might be defunct in a future version
#'   When possible, alternative functions with similar
#'   functionality are also mentioned. 
#' @name pedantics-deprecated
#' @keywords internal
NULL


#' Defunct functions from package \pkg{pedantics}.
#'
#' The functions listed below have not been transferred from pedantics to
#' pedtricks with no short-term plan to implement them since they were relying
#' on C code that need to be heavily updated.
#' 
#' In addition, The questions about rates of pedigree error typically
#' encountered, and their effects on inference of elementary QG parameters are
#' basically settled (Charmantier & Réale 2005, Morrissey *et al.* 2007,
#' Firth *et al.* 2015, Bourret & Garant 2017), and outstanding questions will
#' not plausibly be supported by these functions or any straightforward
#' modifications thereof.
#' 
#' Bourret, A., and D. Garant. 2017. An assessment of the reliability of 
#' quantitative genetics estimates in study systems with high rate of extra
#' pair paternity and low recruitment.  Heredity 118: 229-38.
#'
#' Charmantier, A., and D. Réale. 2005. How do misassigned paternities affect 
#' estimation of heritability in the wild? Molecular Ecology 14: 2839-50
#'
#' Firth, J.A. J.D. Hadfield, A.W. Santure, J. Slate and B.C. Sheldon. 2015.
#' The influence of non-random extra pair paternity on heritability estimates
#'  derived from wild pedigrees. Evolution 69: 1336-44.
#'
#' Morrissey, M.B., A.J. Wilson, J.M. Pemberton, and M.M. Ferguson. 2007. A 
#' framework for power and sensitivity analysis for quantitative genetic studies
#' of natural populations, and a case study in Soay sheep (Ovis aries). Journal
#' of Evolutionary Biology 20: 2309-21.
#'
#' 
#' @name pedantics-defunct
#' @section \code{rpederr}: Permutes a pedigree to create a plausible complete pedigree
#' @section \code{fpederr}: Simulates a pedigree with errors and missing data from a complete pedigree.
#' @keywords internal
NULL

