#' @details
#' * `pedigreeStats()` and `pedStatSummary()` have been replaced by `ped_stats()` with a `summary()` and `plot()` methods to simplify the workflow and allow to get the plots without running the statistics each time
#' * `makePedigreeFactor()` and `makePedigreeNumeric()` have been combined in `convert_ped()` which convert a pedigree to numeric or factor
#' * `fixPedigree()` is now `fix_ped()`
#' * `genomesim()`, `microsim()` and `phensim()` have been renamed as `genome_sim()`, `micro_sim()` and `phen_sim()`
#'
#' For backward compatibility with code using `pedantics` older names are still usable but not recommended.
#'
#' @import dplyr
#' @import tidyr
#' @import ggplot2
#' @import kinship2
#' @importFrom grDevices colors colours dev.new dev.off gray.colors postscript savePlot
#' @importFrom graphics axis barplot hist lines mtext par plot
#' @importFrom stats na.omit rbinom rnorm runif uniroot weighted.mean
#' @importFrom utils flush.console head
#' @importFrom mvtnorm rmvnorm
#' @importFrom grDevices savePlot
#' @importFrom grid grid.segments gpar grid.circle grid.text
#' @importFrom MCMCglmm inverseA rbv
#' @importFrom kinship2 kindepth
#' @importFrom genetics is.genotype nallele allele.names is.locus genotype

"_PACKAGE"

