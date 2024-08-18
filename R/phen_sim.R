#' A function to simulated phenotypic data
#' 
#' @description 
#' Simulates phenotypic data across arbitrary pedigrees.  \code{ phen_sim} simulate direct, maternal and paternal genetic and environmental effects for an arbitrary number of traits with arbitrary patterns of missing data.
#' 
#' @param pedigree A pedigree 
#' @param traits The number of traits for which data should be simulated. 
#' @param randomA An additive genetic covariance matrix, with dimensions a multiple of traits - see details 
#' @param randomE An additive environmental covariance matrix, with dimensions a multiple of traits - see details 
#' @param parentalA A vector indicating which effects in \code{randomA} (if any) to treat as parental effects 
#' @param parentalE A vector indicating which effects in \code{randomE} (if any) to treat as parental effects 
#' @param sampled A vector indicating which individuals are sampled 
#' @param records A single value, array of matrix specifying data record availability - see details 
#' @param returnAllEffects If \code{TRUE} then all individual breeding values and environmental effects are returned 
#' 
#' @details
#' \code{randomA} and \code{randomE} are square matrices with dimension equal to the sum of the number direct and indirect effects.  This must be a multiple of the number of traits, i.e. if an indirect effect is to be simulated for only one of multiple traits, those traits with no indirect effect should be included with (co)variances of zero.
#' 
#' \code{parentalA} and \code{parentalE} are optional vectors of characters indicating which trait positions in \code{randomA} and \code{randomE} are to be treated as indirect effects, and which effects to treat as maternal or paternal.  Valid values are 'd', 'm', and 'p', for direct, maternal indirect and paternal indirect effects, respectively.
#' 
#' \code{records} can be specified either (1) as a single value to be applied to all individuals and traits, (2) as a vector the same length as the number of traits, representing trait-specific rates to be applied uniformly to all individuals, or (3) as data frames with rows for each individual and columns for each trait.  In the third option, observed patterns of data availability can be simulated by supplying 0s and 1s for missing and available individual genotypes, respectively.
#' 
#' 
#' @return
#'   \item{phenotypes}{A dataframe containing phenotypes for all individuals specified to have records.}
#'   \item{allEffects}{(optional) A dataframe with all direct and indirect genetic and environmental effects.}
#' 
#' @seealso \code{\link{micro_sim}}, \code{\link{genome_sim}}
#' 
#' @examples
#' ## make up a pedigree
#' id<-   c("a1","a2","a3","a4","a5","a6","a7","a8","a9")
#' dam<-  c(NA,NA,NA,"a1","a1","a1","a4","a4","a4")
#' sire<- c(NA,NA,NA,"a2","a2","a2","a5","a6","a6")
#' pedigree<-as.data.frame(cbind(id,sire,dam))
#' 
#' traits<-2
#' ## no correlations
#' randomA<-diag(4)
#' randomE<-diag(4)
#' parentalA<-c("d","d","m","m")
#' parentalE<-c("d","d","m","m")
#' 
#' ## generate phenoypic data based on this architecture
#'  phen_sim(pedigree=pedigree,traits=2,randomA=randomA,randomE=randomE,
#'           parentalA=parentalA,parentalE=parentalE)
#' 
#' ## let's do it again but see how the phenotypes were composed
#'  phen_sim(pedigree=pedigree,traits=2,randomA=randomA,randomE=randomE,
#'           parentalA=parentalA,parentalE=parentalE,returnAllEffects=TRUE)
#' 
#'
#' @keywords simulation
#'
#' @export
#'

 phen_sim <-
  function(pedigree, traits = 1, randomA = NULL, randomE = NULL,
           parentalA = NULL, parentalE = NULL, sampled = NULL, records = NULL, returnAllEffects = FALSE) {
    if (is.null(records) == FALSE && is.data.frame(records) == FALSE && is.matrix(records) == FALSE && length(records) != traits && length(records) != 1) {
      stop("Dimension of records incompatible with number of traits")
    }
    if (is.null(records) == FALSE && (is.data.frame(records) | is.matrix(records)) && dim(records)[1] != length(pedigree[, 1])) {
      stop("Dimension of records incompatible with the size of the pedigree")
    }
    if (is.null(records) == FALSE && (is.data.frame(records) | is.matrix(records)) && dim(records)[2] != traits) {
      stop("Dimension of records incompatible with number of traits")
    }


    rmvcn <- function(m, observed, expected, unknown) {
      msize <- dim(m)[1]
      indexes <- as.data.frame(cbind(1:msize, unknown))
      indexes <- indexes[order(indexes[, 2], decreasing = TRUE), ][, 1]

      m <- m[indexes, indexes]
      observed <- observed[indexes]
      expected <- expected[indexes]
      unknowns <- sum(unknown)

      # get the expectation
      u_bar <- expected[1:unknowns] + m[1:unknowns, (unknowns + 1):msize] %*% solve(m[(unknowns + 1):msize, (unknowns + 1):msize]) %*% (observed[(unknowns + 1):msize] - expected[(unknowns + 1):msize])

      # get the Schur complement
      Schur <- m[1:unknowns, 1:unknowns] - m[1:unknowns, (unknowns + 1):msize] %*% solve(m[(unknowns + 1):msize, (unknowns + 1):msize]) %*% m[(unknowns + 1):msize, 1:unknowns]

      # return the expectation and variance
      result <- list(u_bar, Schur, indexes[1:unknowns])
      names(result) <- c("expectation", "variance", "indexes")
      result
    }

    cat(paste("Simulating breeding values and environmental effects..."))
    flush.console()

    n <- length(pedigree[, 1])
    breedingValues <- matrix(0, length(pedigree[, 1]), traits)
    if (is.null(randomA) == FALSE & sum(diag(randomA)) != 0) breedingValues <- rbv(pedigree, randomA)
    environEffects <- NULL
    if (is.matrix(randomE)) {
      environEffects <- rmvnorm(n, rep(0, dim(randomE)[2]), randomE)
    } else {
      environEffects <- rnorm(n, 0, sqrt(randomE))
    }
    effectsTable <- as.data.frame(cbind(pedigree, breedingValues, environEffects))



    for (x in 1:traits) names(effectsTable)[3 + x] <- paste("a_tr", x, sep = "")
    if (is.null(parentalA) == FALSE) {
      for (x in (traits + 1):length(parentalA)) names(effectsTable)[3 + x] <- paste("bv_", parentalA[x], "_tr", x - traits, sep = "")
    }
    if (is.null(parentalA) == FALSE & length(parentalA) == (traits * 3)) {
      for (x in (traits * 2 + 1):length(parentalA)) names(effectsTable)[3 + x] <- paste("bv_", parentalA[x], "_tr", x - traits * 2, sep = "")
    }
    for (x in 1:traits) names(effectsTable)[3 + dim(randomA)[2] + x] <- paste("e_tr", x, sep = "")
    if (is.null(parentalE) == FALSE) {
      for (x in (traits + 1):length(parentalE)) names(effectsTable)[3 + dim(randomA)[2] + x] <- paste("e_", parentalE[x], "_tr", x - traits, sep = "")
    }
    if (is.null(parentalE) == FALSE & length(parentalE) == (traits * 3)) {
      for (x in (traits + 1):length(parentalE)) names(effectsTable)[3 + dim(randomA)[2] + x] <- paste("e_", parentalE[x], "_tr", x - traits * 2, sep = "")
    }


    parentalAeffects <- 0
    if (is.null(parentalA) == FALSE) parentalAeffects <- length(parentalA) - traits
    parentalEeffects <- 0
    if (is.null(parentalE) == FALSE) parentalEeffects <- length(parentalE) - traits

    if (is.null(parentalA) == FALSE) {
      # covariance matrix of parental effects with own effects
      cmmeo <- matrix(NA, 2 * length(parentalA) - traits, 2 * length(parentalA) - traits)
      cmmeo[1:parentalAeffects, 1:parentalAeffects] <- randomA[(traits + 1):length(parentalA), (traits + 1):length(parentalA)]
      cmmeo[(parentalAeffects + 1):(2 * length(parentalA) - traits), (parentalAeffects + 1):(2 * length(parentalA) - traits)] <- randomA
      cmmeo[1:parentalAeffects, (length(parentalA) + 1):(2 * length(parentalA) - traits)] <- 0.25 * randomA[(traits + 1):length(parentalA), (traits + 1):length(parentalA)]
      cmmeo[1:parentalAeffects, (parentalAeffects + 1):length(parentalA)] <- 0.25 * randomA[1:traits, 1:traits]
      cmmeo[(parentalAeffects + 1):(2 * length(parentalA) - traits), 1:parentalAeffects] <- t(cmmeo[1:parentalAeffects, (parentalAeffects + 1):(2 * length(parentalA) - traits)])

      for (x in 1:length(parentalA)) {
        if (parentalA[x] == "m") {
          effectsTable[, (dim(effectsTable)[2] + 1)] <- effectsTable[match(effectsTable$dam, effectsTable[, 1]), (3 + x)]
          names(effectsTable)[dim(effectsTable)[2]] <- paste("P_", names(effectsTable)[(3 + x)], sep = "")
        }
        if (parentalA[x] == "p") {
          effectsTable[, (dim(effectsTable)[2] + 1)] <- effectsTable[match(effectsTable$sire, effectsTable[, 1]), (3 + x)]
          names(effectsTable)[dim(effectsTable)[2]] <- paste("P_", names(effectsTable)[(3 + x)], sep = "")
        }
      }

      for (x in 1:length(pedigree[, 1])) {
        if (is.na(pedigree$dam[x])) {
          obs <- array(dim = 2 * length(parentalA) - traits)
          for (y in 1:length(parentalA)) obs[parentalAeffects + y] <- effectsTable[x, 3 + y]
          exp <- rep(0, (2 * length(parentalA) - traits))
          unk <- c(rep(1, parentalAeffects), rep(0, length(parentalA)))
          parentalDist <- rmvcn(cmmeo, obs, exp, unk)
          parentalvals <- rmvnorm(1, parentalDist$expectation, parentalDist$variance)
          for (y in 1:parentalAeffects) effectsTable[x, (3 + dim(randomA)[2] + dim(randomE)[2] + y)] <- parentalvals[1, y]
        }
      }
    }
    if (is.null(parentalE) == FALSE) {
      for (x in 1:length(parentalE)) {
        if (parentalE[x] == "m") {
          effectsTable[, (dim(effectsTable)[2] + 1)] <- effectsTable[match(effectsTable$dam, effectsTable[, 1]), (3 + dim(randomA)[2] + x)]
          names(effectsTable)[dim(effectsTable)[2]] <- paste("P_", names(effectsTable)[(3 + dim(randomA)[2] + x)], sep = "")
        }
        if (parentalE[x] == "p") {
          effectsTable[, (dim(effectsTable)[2] + 1)] <- effectsTable[match(effectsTable$sire, effectsTable[, 1]), (3 + dim(randomA)[2] + x)]
          names(effectsTable)[dim(effectsTable)[2]] <- paste("P_", names(effectsTable)[(3 + dim(randomA)[2] + x)], sep = "")
        }
      }
      for (x in 1:length(pedigree[, 1])) {
        if (is.na(pedigree$dam[x])) {
          parentalvals <- rmvnorm(1, rep(0, dim(randomE)[2]), randomE)
          for (y in 1:parentalAeffects) effectsTable[x, (3 + dim(randomA)[2] + dim(randomE)[2] + parentalAeffects + y)] <- parentalvals[1, y]
        }
      }
    }

    cat(paste("done.", "\n"))
    flush.console()

    cat(paste("Calculating phenotypes..."))
    flush.console()

    calcPhen <- effectsTable[, 4:length(effectsTable[1, ])]
    retainIndex <- c(rep(1, traits), rep(0, parentalAeffects), rep(1, traits), rep(0, parentalEeffects), rep(1, (parentalAeffects + parentalEeffects)))
    retainIndex <- matrix(retainIndex, length(calcPhen[, 1]), length(calcPhen[1, ]), byrow = TRUE)
    calcPhen <- calcPhen * retainIndex

    for (x in 1:traits) {
      ind <- array(0, dim = traits)
      ind[x] <- 1
      ind <- matrix(ind, length(calcPhen[, 1]), length(calcPhen[1, ]), byrow = TRUE)
      calcPhenTrait <- calcPhen * ind
      effectsTable[, (dim(effectsTable)[2] + 1)] <- rowSums(calcPhenTrait)
      names(effectsTable)[dim(effectsTable)[2]] <- paste("Phen_tr", x, sep = "")
    }

    cat(paste("done.", "\n"))
    flush.console()

    phenAvailIndex <- NULL

    if (is.null(records) == FALSE) {
      cat(paste("Accounting for missing records..."))
      flush.console()
      if (dim(as.data.frame(records))[1] == 1 & dim(as.data.frame(records))[2] == 1) {
        phenAvailIndex <- matrix(rbinom((length(pedigree[, 1]) * traits), 1, records), length(pedigree[, 1]), traits)
      }
      rbinfunc <- function(x) (rbinom(1, 1, x))
      if (traits > 1 & length(records) == traits) {
        recProbs <- matrix(records, length(pedigree[, 1]), traits, byrow = TRUE)
        phenAvailIndex <- apply(recProbs, c(1, 2), rbinfunc)
      }
      if (traits > 1 & length(records) == length(pedigree[, 1])) {
        recProbs <- matrix(records, length(pedigree[, 1]), traits, byrow = FALSE)
        phenAvailIndex <- apply(recProbs, c(1, 2), rbinfunc)
      }
      if (dim(as.data.frame(records))[1] == length(pedigree[, 1]) & dim(as.data.frame(records))[2] == traits) {
        phenAvailIndex <- apply(records, c(1, 2), rbinfunc)
      }
      cat(paste("done.", "\n"))
      flush.console()
    }

    phenotypes <- as.data.frame(effectsTable[, (dim(effectsTable)[2] - traits + 1):dim(effectsTable)[2]])
    phenotypes <- cbind(pedigree[, 1], phenotypes)
    if (is.null(phenAvailIndex) == FALSE) phenotypes <- phenotypes * phenAvailIndex

    names(phenotypes)[1] <- "id"
    for (x in 1:traits) names(phenotypes)[x + 1] <- paste("trait_", x, sep = "")

    if (is.null(sampled) == FALSE) {
      cat(paste("Accounting for unsampled individuals..."))
      flush.console()
      if (is.numeric(sampled) && length(sampled) != length(pedigree[, 1])) {
        stop("Indicator vector for sampled individuals is of different length than pedigree")
      }
      if (is.numeric(sampled)) phenotypes <- subset(phenotypes, sampled > 0)
      if (is.character(sampled)) dataPed <- subset(phenotypes, pedigree[, 1] %in% sampled)
      cat(paste("done.", "\n"))
      flush.console()
    }


    output <- list(phenotypes = phenotypes)
    if (returnAllEffects == TRUE) output <- list(phenotypes = phenotypes, allEffects = effectsTable)

    output
  }

#' @rdname pedantics-deprecated
#' @section \code{phensim}: the function has just been renamed with no other changes for the moment
#' @export
phensim <- function(
  pedigree, traits = 1, randomA = NULL, randomE = NULL,
  parentalA = NULL, parentalE = NULL, sampled = NULL, records = NULL,
  returnAllEffects = FALSE) {
  .Deprecated(genome_sim,
    msg = "this function from pedantics is deprecated, please use the new 'phen_sim()' instead",
  )
  phen_sim(
    pedigree, traits, randomA, randomE,
  parentalA, parentalE, sampled, records,
  returnAllEffects
  )
}
