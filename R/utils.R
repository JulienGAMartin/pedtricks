prune <- function(pedigree, keep, make.base = FALSE) {
  ind.keep <- keep
  nind <- length(ind.keep) + 1
  while (length(ind.keep) != nind) {
    nind <- length(ind.keep)
    ind.keep <- union(na.omit(unlist(pedigree[, 2:3][match(ind.keep, pedigree[, 1]), ])), ind.keep)
  }
  pedigree <- pedigree[sort(match(ind.keep, pedigree[, 1])), ]
  if (make.base) {
    if (any(match(pedigree[, 2], pedigree[, 1]) > match(pedigree[, 1], pedigree[, 1]), na.rm = T)) {
      stop("Dams appearing after their offspring: use fix_ped")
    }
    if (any(match(pedigree[, 3], pedigree[, 1]) > match(pedigree[, 1], pedigree[, 1]), na.rm = T)) {
      stop("Sires appearing after their offspring: use fix_ped")
    }
    phenotyped <- pedigree[, 1] %in% keep
    delete <- rep(FALSE, dim(pedigree)[1])
    for (i in 1:dim(pedigree)[1]) {
      nlinks <- phenotyped[i] + sum(pedigree[, 2] %in% pedigree[, 1][i]) + sum(pedigree[, 3] %in% pedigree[, 1][i]) + sum(is.na(pedigree[i, ][2:3]) == FALSE)
      if (nlinks < 2 & phenotyped[i] == FALSE) {
        pedigree[, 2][which(pedigree[, 2] == pedigree[, 1][i])] <- NA
        pedigree[, 3][which(pedigree[, 3] == pedigree[, 1][i])] <- NA
        delete[i] <- TRUE
      }
    }
    if (any(delete)) {
      pedigree <- pedigree[-which(delete), ]
    }
  }
  pedigree
}

orderPed <- function(ped) {
  reorder <- ped[order(kindepth(ped[, 1], ped[, 2], ped[, 3]), decreasing = FALSE), ]
  return(reorder)
}


getSibNums <- function(ped) {
  ped$damsire <- case_when(
    is.na(ped$dam) ~ NA_character_,
    is.na(ped$sire) ~ NA_character_,
    .default = paste(ped$dam, ped$sire, sep="_")
  )
  f <- table(ped$damsire) * (table(ped$damsire) - 1) / 2 # Number of full siblings
  m <- table(ped$dam) * (table(ped$dam) - 1) / 2 # vector() # Number of maternal siblings
  p <- table(ped$sire) * (table(ped$sire) - 1) / 2 # vector() # Number of paternal siblings
  # Returning results as a vector
  c(full = sum(f), maternal = sum(m), paternal = sum(p))
}

################################################################################
# genotype.list from MasterBayes

#' Genotype Objects for all Loci
#'
#' Creates a \code{list} of \code{genotype} objects from a \code{matrix} or \code{data.frame} of multilocus genotypes.
#'
#' @param G matrix or data.frame of multilocus genotypes with individuals down the rows and loci across columns. Adjacent columns are taken to be the same locus
#' @param marker.type \code{"MSW"} or \code{"MSC"} for co-dominant markers with Wang's (2004) model of genotyping error or CERVUS's model of genotyping error (Kalinowski, 2006; Marshall, 1998) or \code{"AFLP"} for dominant markers (Hadfield, 2009).
#'
#' @return
#'   list of \code{genotype} objects for all loci
#'
#'
#' @author Jarrod Hadfield \email{j.hadfield@ed.ac.uk}
#' @examples
#' \dontrun{
#' data(WarblerG)
#'
#' G <- genotype.list(WarblerG[, -1])
#' summary(G[[1]])
#' }
#'
#' @keywords internal
#' @export

genotype.list <- function(G, marker.type = "MSW") {
  gens <- list()
  if (marker.type == "MSC" | marker.type == "SNP" | marker.type == "MSW") {
    if (length(G[1, ]) %% 2 != 0) {
      stop("Genotypes have odd number of columns")
    }
    for (i in 1:(length(G[1, ]) / 2)) {
      gens[[i]] <- genotype(as.matrix(G[, ((i * 2) - 1):(i * 2)]))
      names(gens)[i] <- names(G[i * 2])
    }
  }
  if (marker.type == "AFLP") {
    for (i in 1:length(G[1, ])) {
      gens[[i]] <- genotypeD(as.matrix(G[, i]))
      names(gens)[i] <- names(G[i])
    }
  }
  gens
}

genotypeD <- function(a1, locus = NULL) {
  alleles <- factor(c("0", "1"))

  object <- factor(as.character(a1), levels = alleles)

  ll <- alleles

  class(object) <- c("genotypeD", "factor")
  attr(object, "allele.names") <- alleles
  attr(object, "allele.map") <- a1

  if (is.null(locus) || is.locus(locus)) {
    attr(object, "locus") <- locus
  } else {
    stop("parameter locus must be of class locus")
  }
  return(object)
}

is.genotypeD <- function(x) {
  inherits(x, "genotypeD")
}

#' @export
print.genotypeD <- function(x, ...) {
  if (!is.null(attr(x, "locus"))) {
    print(attr(x, "locus"))
  }
  print(as.character(x))
  cat("Alleles:", as.character(allele.names(x)), "\n")
  invisible(x)
}

#' @export
summary.genotypeD <- function(object, ...) {
  # if we are called from within summary.data.frame, fall back to
  # summary.factor so that we don't mess up the display


  retval <- list()
  retval$allele.names <- allele.names(object)

  retval$locus <- attr(object, "locus")
  class(retval) <- "summary.genotypeD"

  af <- table(object)

  if (0 %in% names(af)) {
    phat <- af[which(names(af) == 0)] / sum(af[which(names(af) == 0 | names(af) == 1)])
  } else {
    phat <- 0
  }
  p <- sqrt(phat) * (1 / (1 - ((1 - phat) / (8 * sum(af) * phat))))
  q <- 1 - p

  paf <- c(p, q)
  names(paf) <- c("0", "1")

  retval$allele.freq <- cbind("Proportion" = paf)

  gf <- table(object)
  pgf <- prop.table(gf)
  retval$genotype.freq <- cbind("Count" = gf, "Proportion" = pgf)


  ### from code submitted by David Duffy <davidD@qimr.edu.au>
  #
  n.typed <- sum(gf)
  ehet <- (1 - sum(paf * paf))
  matings <- (paf %*% t(paf))^2
  uninf.mating.freq <- sum(matings) - sum(diag(matings))
  pic <- ehet - uninf.mating.freq

  retval$Hu <- 2 * p * q
  retval$pic <- pic
  retval$n.typed <- n.typed
  retval$n.total <- length(object)
  retval$nallele <- nallele(object)
  #
  ###


  if (any(is.na(object))) {
    retval$allele.freq <- rbind(retval$allele.freq, "NA" = c(sum(c(is.na(object), NA))))
    retval$genotype.freq <- rbind(retval$genotype.freq, "NA" = c(sum(c(is.na(object), NA))))
  }

  return(retval)
}

#' @export
`[.genotypeD` <- function(x, i, drop = FALSE) {
  retval <- NextMethod("[")

  # force "NA" not to be a factor level
  ll <- levels(retval) <- na.omit(levels(retval))
  class(retval) <- c("genotypeD", "factor")

  if (drop) {
    alleles <- unique(x)
  } else {
    alleles <- attr(x, "allele.names")
  }

  attr(retval, "allele.names") <- alleles
  attr(retval, "allele.map") <- x
  attr(retval, "locus") <- attr(x, "locus")
  attr(retval, "label") <- attr(x, "label")
  return(retval)
}

#' @export
`[<-.genotypeD` <- function(x, i, value) {
  if (!is.genotypeD(value)) {
    value <- genotypeD(value)
  }

  lx <- levels(x)
  lv <- levels(value)
  ax <- allele.names(x)
  av <- allele.names(value)

  m <- is.na(match(av, ax))
  if (any(m)) {
    warning(paste("Adding new allele name(s):", av[m]))
  }

  la <- unique(c(lx, lv))
  aa <- unique(c(as.character(ax), as.character(av)))

  cx <- class(x)
  nas <- is.na(x)

  data <- match(levels(value)[value], la)

  class(x) <- NULL
  x[i] <- data
  attr(x, "levels") <- la
  map <- attr(x, "allele.map") <- la
  attr(x, "allele.names") <- aa
  class(x) <- cx
  x
}

################################################################################
#'
#' Allele Frequencies
#'
#' Extracts allele frequencies from genotype data
#'
#' @param G data frame or list of \code{genotype} objects
#' @param marker.type \code{"MSW"} or \code{"MSC"} for co-dominant markers with Wang's (2004) model of genotyping error or CERVUS's model of genotyping error (Marshall, 1998) or \code{"AFLP"} for dominant markers.
#'
#' @return list of allele frequencies at each loci
#' @author Jarrod Hadfield \email{j.hadfield@ed.ac.uk}
#' @seealso \code{\link{genotype.list}}, \code{genotype}
#' @examples
#' \dontrun{
#' data(WarblerG)
#'
#' A <- extractA(WarblerG)
#' A[[1]]
#' }
#'
#' @export
#' @keywords internal
#'

extractA <- function(G, marker.type = "MSW") {
  ## from Jarrod to avoid
  ## dependency on MasterBayes

  if (is.genotype(G[[1]]) == FALSE & is.genotypeD(G[[1]]) == FALSE) {
    if ("id" %in% colnames(G)) {
      G <- G[, -which(colnames(G) == "id")]
    }
    if ("categories" %in% colnames(G)) {
      G <- G[, -which(colnames(G) == "categories")]
    }
    G <- genotype.list(G, marker.type = marker.type)
  }

  A <- lapply(G, function(x) {
    summary(x)$allele.freq[, "Proportion"][1:nallele(x)]
  })
  if (marker.type == "AFLP" | is.genotypeD(G[[1]])) {
    A <- lapply(A, function(x) {
      if (any(x == 1)) {
        warning("some loci are monomorphic: setting allele frequencies to 0.9999 and 0.0001")
        x[which(x == 1)] <- 0.9999
        x[which(x == 0)] <- 0.0001
      } else {
        x <- x
      }
      x
    })
  }
  A
}
