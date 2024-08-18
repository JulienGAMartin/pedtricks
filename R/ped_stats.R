#' Calculates a range of statistics of pedigrees
#'
#' Statistics are those that will hopefully be useful for describing pedigrees to be used in quantitative genetic analyses of natural populations.  This module will be most useful when cohort affinities for all individuals can be provided.  All outputs are produced in a numerical form as well as in graphical summaries.
#'
#' @aliases pedigreeStats
#'
#' @param Ped A pedigree
#' @param cohorts (Optional) Cohort affinities for members of the pedigree
#' @param dat (Optional) Available data based upon which the pedigree can be pruned for just informative individuals
#' @param retain The default value ('informative') results in pedigree being pruned to only those individuals who's records contribute to estimation of quantitative genetic parameters with respect to the available data specified in \code{dat}.  Otherwise, specifying a value of 'ancestors' will result in the inclusion of all ancestors of phenotyped individuals.
#' @param includeA If TRUE, additive genetic relatedness matrix is returned.
#' @param lowMem If TRUE, then stats based on calculation of A are not performed.
#'
#' @return
#'   \item{totalMaternities}{Total number of maternities defined by the pedigree.}
#'   \item{totalPaternities}{Total number of paternities defined by the pedigree.}
#'   \item{totalFullSibs}{Total number of pair-wise full sib relationships defined by the pedigree.}
#'   \item{totalMaternalSibs}{Total number of pair-wise maternal sib relationships defined by the pedigree.  To get the number of maternal half sibs, subtract totalFullSibs.}
#'   \item{totalPaternalSibs}{Total number of pair-wise paternal sib relationships defined by the pedigree.  To get the number of paternal half sibs, subtract totalFullSibs.}
#'   \item{totalMaternalGrandmothers}{Total number of maternal grandmothers defined by the pedigree.}
#'   \item{totalMaternalGrandfathers}{Total number of maternal grandfathers defined by the pedigree.}
#'   \item{totalPaternalGrandmothers}{Total number of paternal grandmothers defined by the pedigree.}
#'   \item{totalPaternalGrandfathers}{Total number of paternal grandfathers defined by the pedigree.}
#'   \item{pedigreeDepth}{The pedigree depth, i.e. maximum number of ancestral generations, for each individual.}
#'   \item{inbreedingCoefficients}{Individual inbreeding coefficients}
#'   \item{maternalSibships}{Sibship size of each individual appearing the the dam column of the pedigree.}
#'   \item{paternalSibships}{Sibship size of each individual appearing the the sire column of the pedigree.}
#'   \item{cumulativeRelatedness}{Proportion of pair-wise relatedness values less than values ranging from 0 to 1.}
#'   \item{relatednessCategories}{Discretized distribution of relatedness.}
#'   \item{analyzedPedigree}{Returns the pedigree.}
#'
#'
#'   \item{sampleSizesByCohort}{(Optional) Number of individuals belonging to each cohort.}
#'   \item{maternitiesByCohort}{(Optional) Number of assigned maternities by offspring cohort.}
#'   \item{paternitiesByCohort}{(Optional) Number of assigned paternities by offspring cohort.}
#'   \item{fullSibsByCohort}{(Optional) Number of pair-wise full sib relationships by cohort - note the sum of these need not be equal to totalFullSibs in pedigrees of long-lived organisms.}
#'   \item{maternalSibsByCohort}{(Optional) Number of pair-wise maternal sib relationships by cohort - note the sum of these need not be equal to totalMaternalSibs in pedigrees of long-lived organisms.}
#'   \item{paternalSibsByCohort}{(Optional) Number of pair-wise paternal sib relationships by cohort - note the sum of these need not be equal to totalPaternalSibs in pedigrees of long-lived organisms.}
#'   \item{maternalGrandmothersByCohort}{(Optional) Numbers of maternal grandmother assignments by offspring cohort.}
#'   \item{maternalGrandfathersByCohort}{(Optional) Numbers of maternal grandmother assignments by offspring cohort.}
#'   \item{paternalGrandmothersByCohort}{(Optional) Numbers of paternal grandfather assignments by offspring cohort.}
#'   \item{paternalGrandfathersByCohort}{(Optional) Numbers of paternal grandfather assignments by offspring cohort.}
#'   \item{cumulativePedigreeDepth}{(Optional) Distributions of pedigree depth by cohort.}
#'   \item{meanRelatednessAmongCohorts}{(Optional) Mean relatedness among cohorts.}
#'   \item{cohorts}{(Optional) Returns cohort designations.}
#'
#'   Graphical summaries of a number of these summary statistics are printed to the console when \code{graphicalReports=='y'}.
#'
#' @examples
#' \dontrun{
#'
#' data(gryphons)
#' pedigree <- gryphons[, 1:3]
#'
#' gryphons_ped_stats <- ped_stats(pedigree,
#'   cohorts = gryphons$cohort)
#'
#' gryphons_ped_stats$totalMaternities
#' gryphons_ped_stats$paternitiesByCohort
#' 
#' summary(gryphons_ped_stats)
#'
#' plot(gryphons_ped_stats)

#' }
#' @export
#'


ped_stats <-
  function(Ped, cohorts = NULL, dat = NULL, retain = "informative", includeA = TRUE, lowMem = FALSE) {
    names(Ped)[1] <- "id"
    if (names(Ped)[2] != "dam" | names(Ped)[3] != "sire") {
      if (names(Ped)[3] %in% c("mum", "mom", "mother", "Mum", "Mmom", "Dam", "Mother", "MUM", "MOM", "DAM", "MOTHER")) {
        cat(paste("'mum' appears to be in third column, reordering to 'id','dam','sire'"))
        flush.console()
        Ped <- Ped[, c(1, 3, 2)]
      }
      if (names(Ped)[2] %in% c("mum", "mom", "mother", "Mum", "Mom", "Dam", "Mother", "MUM", "MOM", "DAM", "MOTHER")) {
        names(Ped)[2] <- "dam"
        names(Ped)[3] <- "sire"
      }
      if (names(Ped)[2] != "dam" | names(Ped)[3] != "sire") {
        stop("Unable to identify column names, expecting 'id','dam','sire', or similar")
      }
    }

    for (x in 1:3) Ped[, x] <- as.character(Ped[, x])
    if (is.null(dat) == FALSE && is.numeric(dat) == FALSE) dat <- as.character(dat)

    if (is.null(cohorts) == FALSE && length(Ped[, 1]) != length(cohorts)) {
      stop("Pedigree and cohorts differ in length.")
    }

    if (is.null(dat) == FALSE && is.numeric(dat) && length(Ped[, 1]) != dim(as.data.frame(dat))[1]) {
      stop("Pedigree and available data differ in length.")
    }


    if (is.null(dat) == FALSE) {
      if (is.numeric(dat)) {
        avail <- rowSums(as.data.frame(dat))
        keep <- subset(Ped$id, avail > 0)
      } else {
        keep <- dat
      }

      if (is.null(cohorts) == FALSE) c <- cbind(Ped$id, cohorts)

      if (retain == "informative") {
        Ped <- prune(Ped, keep, make.base = TRUE)
      } else {
        Ped <- prune(Ped, keep, make.base = FALSE)
      }

      if (is.null(cohorts) == FALSE) c <- subset(c, c[, 1] %in% Ped[, 1])
      if (is.null(cohorts) == FALSE) cohorts <- c[, 2]
    }

    # sample size

    totalSampleSize <- length(Ped$id)

    # maternities and paternities

    totalMaternities <- sum(table(Ped$dam))
    totalPaternities <- sum(table(Ped$sire))

    # siblings

    numPed <- convert_ped(type = "numeric", as.character(Ped$id), as.character(Ped$sire),
      as.character(Ped$dam),
      missingVal = NA
    )

    sibNums <- getSibNums(numPed$numericPedigree)

    totalFullSibs <- sibNums["full"]
    totalMaternalSibs <- sibNums["maternal"]
    totalPaternalSibs <- sibNums["paternal"]

    # grandparents

    grandparentData <- Ped
    grandparentData$maternalGM <- grandparentData$dam[match(grandparentData$dam, grandparentData$id)]
    grandparentData$maternalGF <- grandparentData$sire[match(grandparentData$dam, grandparentData$id)]
    grandparentData$paternalGM <- grandparentData$dam[match(grandparentData$sire, grandparentData$id)]
    grandparentData$paternalGF <- grandparentData$sire[match(grandparentData$sire, grandparentData$id)]

    totalMaternalGM <- sum(table(grandparentData$maternalGM))
    totalMaternalGF <- sum(table(grandparentData$maternalGF))
    totalPaternalGM <- sum(table(grandparentData$paternalGM))
    totalPaternalGF <- sum(table(grandparentData$paternalGF))

    # pedigree depth

    pedigreeDepth <- table(kindepth(Ped[, 1], Ped[, 2], Ped[, 3]))

    # inbreeding coefficients

    orderedPed <- as.data.frame(orderPed(Ped))
    orderedPed$inbreeding <- inverseA(orderPed(Ped))$inbreeding
    reorderInbreeding <- as.data.frame(Ped$id)
    reorderInbreeding$inbreeding <- orderedPed$inbreeding[match(reorderInbreeding[, 1], orderedPed$id)]

    # sibship sizes

    matSibships <- as.data.frame(table(as.character(Ped$dam)))
    patSibships <- as.data.frame(table(as.character(Ped$sire)))

    cumulativeRelatedness <- NULL
    pairwiseRelatedness <- NULL
    relatednessBin <- NULL
    if (lowMem == FALSE) {
      # relatedness classes
      cutoffs <- seq(-0.0125, 0.9875, by = 0.025)
      midBins <- seq(0, 0.975, by = 0.025)
      cumulativeRelatedness <- array(dim = length(midBins))
      names(cumulativeRelatedness) <- midBins
      relatednessBin <- array(dim = length(midBins))
      names(relatednessBin) <- midBins
      A <- kinship(Ped[, 1], Ped[, 3], Ped[, 2]) * 2
      pairwiseRelatedness <- A * abs(diag(length(A[, 1])) - 1)
      pairwiseRelatedness[upper.tri(pairwiseRelatedness)] <- 0
      pairwiseRelatedness <- c(pairwiseRelatedness)
      #  pairwiseRelatedness<-subset(pairwiseRelatedness,pairwiseRelatedness>1e-9)
      for (x in 2:(length(cutoffs) - 1)) {
        #  cumulativeRelatedness[x]<-table(pairwiseRelatedness<cutoffs[x+1])["TRUE"]
        relatednessBin[x] <- table(pairwiseRelatedness > cutoffs[x] & pairwiseRelatedness <= cutoffs[x + 1])["TRUE"]
      }
      relatednessBin[1] <- ((totalSampleSize^2 - totalSampleSize) / 2) - sum(relatednessBin, na.rm = TRUE)
      relatednessBin[relatednessBin %in% NA] <- 0

      rb <- relatednessBin / sum(relatednessBin)
      for (x in 1:(length(cutoffs) - 1)) {
        cumulativeRelatedness[x] <- sum(rb[1:x])
      }
    }

    # MacCluer's pedigree completeness statistics
    missingness <- NULL
    missingness <- as.data.frame(cbind(as.character(Ped$id), matrix(0, length(Ped[, 1]), max(kindepth(Ped[, 1], Ped[, 2], Ped[, 3])))))
    p <- Ped
    for (x in 1:(length(missingness[1, ]) - 1)) {
      pOld <- p
      if (is.null(pOld[, 2]) == FALSE & is.null(pOld[, 3]) == FALSE) {
        matsub <- subset(p, is.na(p[, 2]) == FALSE)[, c(1, 2)]
        names(matsub) <- c("id", "p")
        patsub <- subset(p, is.na(p[, 3]) == FALSE)[, c(1, 3)]
        names(patsub) <- c("id", "p")
        p <- as.data.frame(rbind(matsub, patsub))
      }
      if (is.null(pOld[, 2]) == FALSE & is.null(pOld[, 3]) == TRUE) {
        p <- as.data.frame(subset(p, is.na(p[, 2]) == FALSE)[, c(1, 2)])
      }
      if (is.null(pOld[, 2]) == TRUE & is.null(pOld[, 3]) == FALSE) {
        p <- as.data.frame(subset(p, is.na(p[, 3]) == FALSE)[, c(1, 3)])
      }
      t <- table(p[, 1])
      missingness[, 1 + x] <- t[match(missingness[, 1], names(t))]
      p$gm <- Ped[match(p[, 2], Ped[, 1]), 2]
      p$gf <- Ped[match(p[, 2], Ped[, 1]), 3]
      p <- p[, c(1, 3, 4)]
    }
    missingness[is.na(missingness)] <- 0
    for (x in 1:(length(missingness[1, ]) - 1)) missingness[, 1 + x] <- missingness[, 1 + x] / (2^x)


    ##### if cohort designations are available

    if (is.null(cohorts) == FALSE) {
      # average relatedness within and between cohorts
      vector.to.design.matrix <- function(f) {
        factors <- as.numeric(as.factor(f))
        design.matrix <- matrix(0, length(factors), nlevels(as.factor(f)))
        for (x in 1:length(factors)) design.matrix[x, factors[x]] <- 1
        results <- list(design.matrix, levels(as.factor(f)))
        names(results) <- c("DesignMatrix", "FactorLevels")
        results
      }
      Y <- vector.to.design.matrix(cohorts)$DesignMatrix
      D <- diag(1 / table(cohorts))

      meanRelatednessAmongCohorts <- NULL
      if (lowMem == FALSE) {
        meanRelatednessAmongCohorts <- t(Y %*% D) %*% A %*% (Y %*% D)
        rownames(meanRelatednessAmongCohorts) <- names(table(cohorts))
        colnames(meanRelatednessAmongCohorts) <- names(table(cohorts))
      }

      # sample size by cohort

      cohortSampleSizes <- table(cohorts)

      # maternities and paternities by cohort

      cohortMaternities <- array(0, dim = length(table(cohorts)))
      cohortPaternities <- array(0, dim = length(table(cohorts)))
      names(cohortMaternities) <- names(table(cohorts))
      names(cohortPaternities) <- names(table(cohorts))

      cohortFullSibs <- array(0, dim = length(table(cohorts)))
      cohortMaternalSibs <- array(0, dim = length(table(cohorts)))
      cohortPaternalSibs <- array(0, dim = length(table(cohorts)))

      cohortMaternalGM <- array(0, dim = length(table(cohorts)))
      cohortMaternalGF <- array(0, dim = length(table(cohorts)))
      cohortPaternalGM <- array(0, dim = length(table(cohorts)))
      cohortPaternalGF <- array(0, dim = length(table(cohorts)))
      names(cohortMaternalGM) <- names(table(cohorts))
      names(cohortMaternalGF) <- names(table(cohorts))
      names(cohortPaternalGM) <- names(table(cohorts))
      names(cohortPaternalGF) <- names(table(cohorts))

      cohortPedgireeDepth <- matrix(0, length(table(cohorts)), max(as.numeric(names(pedigreeDepth))) + 1)
      rownames(cohortPedgireeDepth) <- names(table(cohorts))
      colnames(cohortPedgireeDepth) <- as.character(0:max(as.numeric(names(pedigreeDepth))))

      for (x in 1:length(cohortMaternities)) {
        temp <- subset(Ped, as.character(cohorts) == names(table(cohorts))[x])
        cohortMaternities[x] <- sum(table(temp$dam))
        cohortPaternities[x] <- sum(table(temp$sire))

        numPed <- convert_ped(type = "numeric", as.character(temp$id), as.character(temp$sire),
          as.character(temp$dam),
          missingVal = NA
        )
        sibNums <- getSibNums(numPed$numericPedigree)
        cohortFullSibs[x] <- sibNums["full"]
        cohortMaternalSibs[x] <- sibNums["maternal"]
        cohortPaternalSibs[x] <- sibNums["paternal"]

        temp <- subset(Ped, as.numeric(as.character(cohorts)) <= as.numeric(names(table(cohorts))[x]))
        pedDepth <- table(kindepth(temp[, 1], temp[, 2], temp[, 3]))
        for (y in 1:length(pedDepth)) cohortPedgireeDepth[x, names(pedDepth[y])] <- pedDepth[y]

        temp <- subset(grandparentData, as.character(cohorts) == names(table(cohorts))[x])
        cohortMaternalGM[x] <- sum(table(temp$maternalGM))
        cohortMaternalGF[x] <- sum(table(temp$maternalGF))
        cohortPaternalGM[x] <- sum(table(temp$paternalGM))
        cohortPaternalGF[x] <- sum(table(temp$paternalGF))
      }
      cohortFullSibs <- as.data.frame(cohortFullSibs)
      cohortMaternalSibs <- as.data.frame(cohortMaternalSibs)
      cohortPaternalSibs <- as.data.frame(cohortPaternalSibs)
      rownames(cohortFullSibs) <- names(table(cohorts))
      rownames(cohortMaternalSibs) <- names(table(cohorts))
      rownames(cohortPaternalSibs) <- names(table(cohorts))


      matSibships$mumBirthYear <- cohorts[match(matSibships[, 1], Ped$id)]
      patSibships$mumBirthYear <- cohorts[match(patSibships[, 1], Ped$id)]



      # MacCluer's pedigree completeness statistics, averaged by cohort
      cohortPedigreeCompleteness <- NULL
      cohortPedigreeCompleteness <- matrix(NA, length(table(cohorts)), length(missingness[1, ]) - 1)
      for (x in 1:length(table(cohorts))) {
        subsetMissingness <- subset(missingness, as.character(cohorts) == names(table(cohorts))[x])
        for (y in 1:length(cohortPedigreeCompleteness[1, ])) {
          cohortPedigreeCompleteness[x, y] <- mean(subsetMissingness[, y + 1])
        }
      }
    }

    if (lowMem == FALSE & includeA == TRUE) {
      results <- list(
        totalSampleSize = totalSampleSize,
        totalMaternities = totalMaternities,
        totalPaternities = totalPaternities,
        totalFullSibs = totalFullSibs[[1]],
        totalMaternalSibs = totalMaternalSibs[[1]],
        totalPaternalSibs = totalPaternalSibs[[1]],
        totalMaternalGrandmothers = totalMaternalGM,
        totalMaternalGrandfathers = totalMaternalGF,
        totalPaternalGrandmothers = totalPaternalGM,
        totalPaternalGrandfathers = totalPaternalGF,
        pedigreeDepth = pedigreeDepth,
        inbreedingCoefficients = reorderInbreeding$inbreeding,
        Amatrix = A,
        maternalSibships = matSibships,
        paternalSibships = patSibships,
        cumulativeRelatedness = cumulativeRelatedness,
        relatednessCategories = relatednessBin,
        analyzedPedigree = Ped,
        missingness = missingness,
        cohortPedigreeCompleteness = cohortPedigreeCompleteness
      )

      if (is.null(cohorts) == FALSE) {
        results <- list(
          totalSampleSize = totalSampleSize,
          totalMaternities = totalMaternities,
          totalPaternities = totalPaternities,
          totalFullSibs = totalFullSibs[[1]],
          totalMaternalSibs = totalMaternalSibs[[1]],
          totalPaternalSibs = totalPaternalSibs[[1]],
          totalMaternalGrandmothers = totalMaternalGM,
          totalMaternalGrandfathers = totalMaternalGF,
          totalPaternalGrandmothers = totalPaternalGM,
          totalPaternalGrandfathers = totalPaternalGF,
          sampleSizesByCohort = cohortSampleSizes,
          maternitiesByCohort = cohortMaternities,
          paternitiesByCohort = cohortPaternities,
          fullSibsByCohort = cohortFullSibs,
          maternalSibsByCohort = cohortMaternalSibs,
          paternalSibsByCohort = cohortPaternalSibs,
          maternalGrandmothersByCohort = cohortMaternalGM,
          maternalGrandfathersByCohort = cohortMaternalGF,
          paternalGrandmothersByCohort = cohortPaternalGM,
          paternalGrandfathersByCohort = cohortPaternalGF,
          pedigreeDepth = pedigreeDepth,
          cumulativePedigreeDepth = cohortPedgireeDepth,
          meanRelatednessAmongCohorts = meanRelatednessAmongCohorts,
          inbreedingCoefficients = reorderInbreeding$inbreeding,
          Amatrix = A,
          maternalSibships = matSibships,
          paternalSibships = patSibships,
          cumulativeRelatedness = cumulativeRelatedness,
          relatednessCategories = relatednessBin,
          cohorts = cohorts,
          analyzedPedigree = Ped,
          missingness = missingness,
          cohortPedigreeCompleteness = cohortPedigreeCompleteness
        )
      }
    } else {
      results <- list(
        totalSampleSize = totalSampleSize,
        totalMaternities = totalMaternities,
        totalPaternities = totalPaternities,
        totalFullSibs = totalFullSibs[[1]],
        totalMaternalSibs = totalMaternalSibs[[1]],
        totalPaternalSibs = totalPaternalSibs[[1]],
        totalMaternalGrandmothers = totalMaternalGM,
        totalMaternalGrandfathers = totalMaternalGF,
        totalPaternalGrandmothers = totalPaternalGM,
        totalPaternalGrandfathers = totalPaternalGF,
        pedigreeDepth = pedigreeDepth,
        inbreedingCoefficients = reorderInbreeding$inbreeding,
        maternalSibships = matSibships,
        paternalSibships = patSibships,
        cumulativeRelatedness = cumulativeRelatedness,
        relatednessCategories = relatednessBin,
        analyzedPedigree = Ped,
        missingness = missingness
      )

      if (is.null(cohorts) == FALSE) {
        results <- list(
          totalSampleSize = totalSampleSize,
          totalMaternities = totalMaternities,
          totalPaternities = totalPaternities,
          totalFullSibs = totalFullSibs[[1]],
          totalMaternalSibs = totalMaternalSibs[[1]],
          totalPaternalSibs = totalPaternalSibs[[1]],
          totalMaternalGrandmothers = totalMaternalGM,
          totalMaternalGrandfathers = totalMaternalGF,
          totalPaternalGrandmothers = totalPaternalGM,
          totalPaternalGrandfathers = totalPaternalGF,
          sampleSizesByCohort = cohortSampleSizes,
          maternitiesByCohort = cohortMaternities,
          paternitiesByCohort = cohortPaternities,
          fullSibsByCohort = cohortFullSibs,
          maternalSibsByCohort = cohortMaternalSibs,
          paternalSibsByCohort = cohortPaternalSibs,
          maternalGrandmothersByCohort = cohortMaternalGM,
          maternalGrandfathersByCohort = cohortMaternalGF,
          paternalGrandmothersByCohort = cohortPaternalGM,
          paternalGrandfathersByCohort = cohortPaternalGF,
          pedigreeDepth = pedigreeDepth,
          cumulativePedigreeDepth = cohortPedgireeDepth,
          meanRelatednessAmongCohorts = meanRelatednessAmongCohorts,
          inbreedingCoefficients = reorderInbreeding$inbreeding,
          maternalSibships = matSibships,
          paternalSibships = patSibships,
          cumulativeRelatedness = cumulativeRelatedness,
          relatednessCategories = relatednessBin,
          cohorts = cohorts,
          analyzedPedigree = Ped,
          missingness = missingness,
          cohortPedigreeCompleteness = cohortPedigreeCompleteness
        )
      }
    }
    class(results) <- "ped_stats"
    results
  }

#' @rdname pedantics-deprecated
#' @section \code{pedigreeStats}: the function has been simplified but only the functionality are still available via \code{ped_stats} and its summary and plot methods
#' @export
pedigreeStat <- function() {
  .Deprecated(ped_stats,
    msg = "this function from pedantics is deprecated and not working anymore. Please use 'ped_stats()' instead",
  )
}
