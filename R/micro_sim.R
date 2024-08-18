#' Simulates microsatellite data across a pedigree.
#'
#' Uses a pedigree with parents identified for all non-founding
#' individuals and simulates microsatellite genotypes
#'
#' @param pedigree A pedigree
#' @param genFreqs (optional) A list of allele frequencies, can be produced with \code{extractA}
#' @param genotypesSample (required if \code{genFreqs} is not supplied) a sample of genotypes from which to estimate population allele frequencies
#' @param knownGenotypes (not yet implemented) a data frame of genotypes for (potentially a subset) of founder individuals
#' @param records Record availability, see details.
#' @param eRate1 The rate of genotypic substitution errors, i.e., when a true genotype at a given locus is replaced by a pair of alleles selected at random based on the population allele frequencies
#' @param eRate2 The rate of allelic substitution errors, i.e. when an allele is erroneously replaced at a given locus by an allele chosen at random based on the population allele frequencies
#' @param eRate3 The rate of large allele dropouts, simulated by setting the value of the larger allele at a locus to the value of the smaller allele
#'
#' @details
#'   Error rates and data availability rates can be specified as either (1) single values to be applied to all individuals and all loci, (2) as a vector the same length as the number of loci, representing locus-specific rates to be applied uniformly to all individuals, or (3) as data frames with rows for each individual and columns for each locus.  In the third option, observed patterns of data availability can be simulated by supplying 0s and 1s for missing and available individual genotypes, respectively.
#'
#' @return
#'   \item{trueGenotypes}{A data frame of true genotypes}
#'   \item{observedGenotypes}{A data frame of plausible observed genotypes, given specified patterns of missingness and errors.}
#'
#'
#' @seealso \code{\link{phen_sim}}, \code{\link{genome_sim}}
#'
#' @examples
#'
#' pedigree <- as.data.frame(matrix(c(
#'   "m1",   NA,     NA,
#'   "m2",   NA,     NA,
#'   "m3",   NA,     NA,
#'   "d4",   NA,     NA,
#'   "d5",   NA,     NA,
#'   "o6",   "m1",   "d4",
#'   "o7",   "m1",   "d4",
#'   "o8",   "m1",   "d4",
#'   "o9",   "m1",   "d4",
#'   "o10",  "m2",   "d5",
#'   "o11",  "m2",   "d5",
#'   "o12",  "m2",   "d5",
#'   "o13",  "m2",   "d5",
#'   "o14",  "m3",   "d5",
#'   "o15",  "m3",   "d5",
#'   "o16",  "m3",   "d5",
#'   "o17",  "m3",   "d5"
#' ), 17, 3, byrow = TRUE))
#' names(pedigree) <- c("id", "dam", "sire")
#' for (x in 1:3) pedigree[, x] <- as.factor(pedigree[, x])
#'
#' ## some sample genotypes, very simple, two markers with He = 0.5
#' sampleGenotypes <- as.data.frame(matrix(c(
#'   1, 2, 1, 2, 2, 1, 2, 1
#' ), 2, 4, byrow = TRUE))
#' ## locus names
#' names(sampleGenotypes) <- c("loc1a", "loc1b", "loc2a", "loc2b")
#'
#' ## simulate some genotypes
#' micro_sim(pedigree = pedigree, genotypesSample = sampleGenotypes)
#'
#' @keywords simulation
#'
#' @export
#'



micro_sim <-
  function(pedigree, genFreqs = NULL, genotypesSample = NULL, knownGenotypes = NULL, records = NULL, eRate1 = 0, eRate2 = 0, eRate3 = 0) {
    if (is.null(genotypesSample) == FALSE && (names(genotypesSample)[1] == "id" | names(genotypesSample)[1] == "ID")) genotypesSample <- genotypesSample[, -1]
    if (is.null(genotypesSample) == FALSE & is.null(genFreqs) == TRUE) {
      genFreqs <- extractA(genotypesSample)
    }

    loci <- length(names(genFreqs))

    if (is.null(records) == FALSE && is.data.frame(records) == FALSE && is.matrix(records) == FALSE && length(records) != loci && length(records) != 1) {
      stop("Dimension of records incompatible with number of loci")
    }
    if (is.null(records) == FALSE && (is.data.frame(records) | is.matrix(records)) && dim(records)[1] != length(pedigree[, 1])) {
      stop("Dimension of records incompatible with the size of the pedigree")
    }
    if (is.null(records) == FALSE && (is.data.frame(records) | is.matrix(records)) && dim(records)[2] != loci) {
      stop("Dimension of records incompatible with number of loci")
    }

    if (length(eRate1) != 1 & eRate1 != loci) stop("Length of eRate1 incompatible with number of loci")
    if (length(eRate2) != 1 & eRate2 != loci) stop("Length of eRate2 incompatible with number of loci")
    if (length(eRate3) != 1 & eRate3 != loci) stop("Length of eRate3 incompatible with number of loci")

    microGen <- as.data.frame(matrix(NA, length(pedigree[, 1]), loci * 2 + 1))
    microGen[, 1] <- as.character(pedigree[, 1])
    row.names(microGen) <- as.character(pedigree[, 1])
    microGen$genotyped <- 0

    randAllele <- function(loc, freqs) {
      i <- runif(1, 0, 1)
      a <- 0
      while (i > 0) {
        a <- a + 1
        i <- i - freqs[[loc]][a]
      }
      names(freqs[[loc]][a])
    }


    assignGenotype <- function(id, mum, dad, f) {
      maternalGenotype <- NULL
      paternalGenotype <- NULL

      if (is.na(mum) == TRUE) {
        maternalGenotype <- array(dim = loci * 2)
        for (l in 1:loci) {
          maternalGenotype[(l - 1) * 2 + 1] <- randAllele(l, f)
          maternalGenotype[(l - 1) * 2 + 2] <- randAllele(l, f)
        }
      } else {
        if (microGen[mum, dim(microGen)[2]] == 1) {
          maternalGenotype <- subset(microGen, microGen[, 1] == mum)[1, 2:(loci * 2 + 1)]
        }
      }

      if (is.na(dad) == TRUE) {
        paternalGenotype <- array(dim = loci * 2)
        for (l in 1:loci) {
          paternalGenotype[(l - 1) * 2 + 1] <- randAllele(l, f)
          paternalGenotype[(l - 1) * 2 + 2] <- randAllele(l, f)
        }
      } else {
        if (microGen[dad, dim(microGen)[2]] == 1) {
          paternalGenotype <- subset(microGen, microGen[, 1] == dad)[1, 2:(loci * 2 + 1)]
        }
      }

      if (is.null(maternalGenotype) == FALSE & is.null(paternalGenotype) == FALSE) {
        microGen[id, dim(microGen)[2]] <<- 1
        for (l in 1:loci) {
          maternalAllele <- rbinom(1, 1, 0.5)
          paternalAllele <- rbinom(1, 1, 0.5)
          microGen[id, (l - 1) * 2 + 2] <<- maternalGenotype[(l - 1) * 2 + 1 + maternalAllele]
          microGen[id, (l - 1) * 2 + 3] <<- paternalGenotype[(l - 1) * 2 + 1 + paternalAllele]
        }
      }
    }

    while (prod(microGen$genotyped) == 0) {
      for (x in 1:dim(microGen)[1]) {
        assignGenotype(as.character(pedigree$id[x]), as.character(pedigree$dam[x]), as.character(pedigree$sire[x]), genFreqs)
      }
    }

    trueMicroGenotypes <- microGen[-c(1, dim(microGen)[2])]
    for (l in 1:loci) names(trueMicroGenotypes)[(l - 1) * 2 + 1] <- names(genFreqs[l])
    for (l in 1:loci) names(trueMicroGenotypes)[(l - 1) * 2 + 2] <- paste(names(genFreqs[l]), "_b", sep = "")


    observedMicroGenotypes <- trueMicroGenotypes

    # make some genotypic errors

    # genotype substitutions
    if (sum(eRate1) > 0) {
      if (length(eRate1) == 1) {
        for (i in 1:dim(observedMicroGenotypes)[1]) {
          for (l in 1:loci) {
            if (runif(1, 0, 1) < eRate1) {
              observedMicroGenotypes[i, (l - 1) * 2 + 1] <- randAllele(l, genFreqs)
              observedMicroGenotypes[i, (l - 1) * 2 + 2] <- randAllele(l, genFreqs)
            }
          }
        }
      }
      if (length(eRate1) == loci) {
        for (i in 1:dim(observedMicroGenotypes)[1]) {
          for (l in 1:loci) {
            if (runif(1, 0, 1) < eRate1[l]) {
              observedMicroGenotypes[i, (l - 1) * 2 + 1] <- randAllele(l, genFreqs)
              observedMicroGenotypes[i, (l - 1) * 2 + 2] <- randAllele(l, genFreqs)
            }
          }
        }
      }
    }
    # allele substitutions
    if (sum(eRate2) > 0) {
      if (length(eRate2) == 1) {
        for (i in 1:dim(observedMicroGenotypes)[1]) {
          for (l in 1:loci) {
            if (runif(1, 0, 1) < eRate2) observedMicroGenotypes[i, (l - 1) * 2 + 1] <- randAllele(l, genFreqs)
            if (runif(1, 0, 1) < eRate2) observedMicroGenotypes[i, (l - 1) * 2 + 2] <- randAllele(l, genFreqs)
          }
        }
      }
      if (length(eRate2) == loci) {
        for (i in 1:dim(observedMicroGenotypes)[1]) {
          for (l in 1:loci) {
            if (runif(1, 0, 1) < eRate2[l]) observedMicroGenotypes[i, (l - 1) * 2 + 1] <- randAllele(l, genFreqs)
            if (runif(1, 0, 1) < eRate2[l]) observedMicroGenotypes[i, (l - 1) * 2 + 2] <- randAllele(l, genFreqs)
          }
        }
      }
    }
    # large allele dropouts
    if (sum(eRate3) > 0) {
      if (length(eRate3) == 1) {
        for (i in 1:dim(observedMicroGenotypes)[1]) {
          for (l in 1:loci) {
            if (runif(1, 0, 1) < eRate3) {
              smallerAllele <- min(c(trueMicroGenotypes[i, (l - 1) * 2 + 1], trueMicroGenotypes[i, (l - 1) * 2 + 2]))
              observedMicroGenotypes[i, (l - 1) * 2 + 1] <- smallerAllele
              observedMicroGenotypes[i, (l - 1) * 2 + 2] <- smallerAllele
            }
          }
        }
      }
      if (length(eRate3) == loci) {
        for (i in 1:dim(observedMicroGenotypes)[1]) {
          for (l in 1:loci) {
            if (runif(1, 0, 1) < eRate3[l]) {
              smallerAllele <- min(c(trueMicroGenotypes[i, (l - 1) * 2 + 1], trueMicroGenotypes[i, (l - 1) * 2 + 2]))
              observedMicroGenotypes[i, (l - 1) * 2 + 1] <- smallerAllele
              observedMicroGenotypes[i, (l - 1) * 2 + 2] <- smallerAllele
            }
          }
        }
      }
    }

    # cut out some missing records
    genAvailIndex <- NULL
    if (is.null(records) == FALSE) {
      recProbs <- matrix(0, dim(trueMicroGenotypes)[1], dim(trueMicroGenotypes)[2])
      if (dim(as.data.frame(records))[1] == 1 & dim(as.data.frame(records))[2] == 1) {
        recProbs <- matrix(records, dim(trueMicroGenotypes)[1], dim(trueMicroGenotypes)[2])
      }
      if (length(records) == loci) for (l in 1:loci) recProbs[, (l - 1) * 2 + 1] <- records[l]
      if (dim(as.data.frame(records))[1] == length(pedigree[, 1]) & dim(as.data.frame(records))[2] == loci) {
        for (l in 1:loci) recProbs[, (l - 1) * 2 + 1] <- records[, l]
      }
      rbinfunc <- function(x) (rbinom(1, 1, x))
      genAvailIndex <- apply(recProbs, c(1, 2), rbinfunc)
      for (l in 1:loci) genAvailIndex[, (l - 1) * 2 + 2] <- genAvailIndex[, (l - 1) * 2 + 1]
    }

    if (is.null(genAvailIndex) == FALSE) {
      missingtozeros <- matrix(as.numeric(as.matrix(observedMicroGenotypes)), dim(observedMicroGenotypes)[1], dim(observedMicroGenotypes)[2], byrow = FALSE) * genAvailIndex
      zeromissingfunc <- function(x) {
        res <- x
        if (x == 0) res <- NA
        res
      }
      withMissing <- as.data.frame(apply(missingtozeros, c(1, 2), zeromissingfunc))
      names(withMissing) <- names(observedMicroGenotypes)
      row.names(withMissing) <- row.names(observedMicroGenotypes)
      observedMicroGenotypes <- withMissing
    }

    list(trueGenotypes = trueMicroGenotypes, observedGenotypes = observedMicroGenotypes)
  }

#' @rdname pedantics-deprecated
#' @section \code{microsim}: the function has just been renamed with no other changes for the moment
#' @export
microsim <- function(pedigree, genFreqs = NULL, genotypesSample = NULL, knownGenotypes = NULL, records = NULL, eRate1 = 0, eRate2 = 0, eRate3 = 0) {
  .Deprecated(genome_sim,
    msg = "this function from pedantics is deprecated, please use the new 'micro_sim()' instead",
  )
  micro_sim(
    pedigree, genFreqs, genotypesSample, knownGenotypes, records, eRate1, eRate2, eRate3
  )
}
