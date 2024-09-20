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


#####################
### Cousin Functions
#####################

## function to make products of all combinations and sum them 
combo_prod <- function(x) if(length(x)>1){sum( utils::combn(x, m =2)[1, ] * utils::combn(x, m =2)[2, ])}else{0}

## function to get cousins via a certain parent (p) and grandparent (gp) from a certain pedigree (ped)
## assumes that you have put grandparents into the pedigree  and grandparent pairs (so 9 columns) 
n_cousin <- function(p,gp,ped){
  if(all(is.na(ped[,gp]))){ ## stops working when there no links through a certain grandparent type
    0
  }else{
    f1 <- formula(paste("id~",p,"+",gp))
    f2 <- formula(paste("id~",gp))
    d1 <- aggregate(f1, ped,length)
    sum(aggregate(f2,d1,combo_prod)$id) 
  }
}

## function to work out all the cousin relationships
getCousinNums <- function(ped) {

  cousin_D_FS <- n_cousin(p="dam", gp="maternalGP", ped=ped)
  cousin_D_MHS <- n_cousin(p="dam", gp="maternalGM", ped=ped) - cousin_D_FS
  cousin_D_PHS <- n_cousin(p="dam", gp="maternalGF", ped=ped) - cousin_D_FS

  cousin_S_FS <- n_cousin(p="sire", gp="paternalGP", ped=ped)
  cousin_S_MHS <- n_cousin(p="sire", gp="paternalGM", ped=ped) - cousin_S_FS
  cousin_S_PHS <- n_cousin(p="sire", gp="paternalGF", ped=ped) - cousin_S_FS


    ## stack maternal and paternal grandparents, to get allows for all cousins relationships, not just through mothers or father
  mat_gp <- ped[,c("id","dam","maternalGM","maternalGF","maternalGP")]
  pat_gp <- ped[,c("id","sire","paternalGM","paternalGF","paternalGP")]
  colnames(pat_gp) <- colnames(mat_gp)

  GM <- rbind(pat_gp,mat_gp)
  
  cousin_FS <- n_cousin(p="dam", gp="maternalGP", ped=GM)
  cousin_DS_FS <- cousin_FS-cousin_D_FS-cousin_S_FS

  cousin_MS <- n_cousin(p="dam", gp="maternalGM", ped=GM)
  cousin_DS_MHS <- cousin_MS-cousin_FS-cousin_D_MHS-cousin_S_MHS

  cousin_PS <- n_cousin(p="dam", gp="maternalGF", ped=GM)
  cousin_DS_PHS <- cousin_PS-cousin_FS-cousin_D_PHS-cousin_S_PHS
  
## double cousins - share both sets of gps

  c(cousin_D_FS=cousin_D_FS,
    cousin_DS_FS=cousin_DS_FS,
    cousin_S_FS=cousin_S_FS,
    cousin_D_MHS=cousin_D_MHS,
    cousin_DS_MHS=cousin_DS_MHS,
    cousin_S_MHS=cousin_S_MHS,
    cousin_D_PHS=cousin_D_PHS,
    cousin_DS_PHS=cousin_DS_PHS,
    cousin_S_PHS=cousin_S_PHS)
# c(cousin_D_FS=cousin_D_FS,
#     cousin_DS_FS=cousin_DS_FS,
#     cousin_S_FS=cousin_S_FS,
#     cousin_D_HS=cousin_D_MHS + cousin_D_PHS,
#     cousin_DS_HS=cousin_DS_MHS + cousin_DS_PHS,
#     cousin_S_HS=cousin_S_MHS + cousin_S_PHS)



}
#####################
### Aunt/uncle Functions
#####################

## function to get aunts and uncles via a certain  grandparent (link) from a certain pedigree (ped)
## assumes that you have put grandparents and grandparent pairs into the pedigree (so 9 columns) 
n_au <- function(link,ped){
  if(all(is.na(ped[,link]))){
    0
  }else{  
    p <- if(grepl("maternal",link)){"dam"}else{"sire"}
    gp <- if(grepl("M",link)){"dam"}else if(grepl("F",link)){"sire"} else{"pair"}
    f1 <- formula(paste("id~",p,"+",link))
    f2 <- formula(paste("id~",gp))
    d1 <- aggregate(f1, ped,length)
    d2 <- aggregate(f2, ped,length)
    sum(d1$id * (d2[match(d1[,link],d2[,gp]),"id"]-1), na.rm=TRUE)
  }
}

## function to work out all the aunt/uncle relationships

getAuNums <- function(ped){

  # related through a maternal grandmother (au_D_FS + au_D_MHS ?)
  au_D_MS <- n_au("maternalGM",ped)

  # related through a maternal grandfather (au_D_FS + au_D_PHS ?)
  au_D_PS <- n_au("maternalGF",ped)

  # related through a maternal grandparents (au_D_FS ?)
  au_D_FS <- n_au("maternalGP",ped)

  au_D_PHS <- au_D_PS - au_D_FS
  au_D_MHS <- au_D_MS - au_D_FS

  # related through a paternal grandmother (au_S_FS + au_S_MHS ?)
  au_S_MS <- n_au("paternalGM",ped)

  # related through a paternal grandfather (au_S_FS + au_D_PHS ?)
  au_S_PS <- n_au("paternalGF",ped)

  # related through a paternal grandparents (au_S_FS ?)
  au_S_FS <- n_au("paternalGP",ped)

  au_S_PHS <- au_S_PS - au_S_FS
  au_S_MHS <- au_S_MS - au_S_FS

  c(au_D_FS=au_D_FS,au_S_FS=au_S_FS,au_D_MHS=au_D_MHS,au_S_MHS=au_S_MHS,au_D_PHS=au_D_PHS,au_S_PHS=au_S_PHS)
}



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

### exported function providing methods

#' @export
print.genotypeD <- function(x, ...) {
  if (!is.null(attr(x, "locus"))) {
    print(attr(x, "locus"))
  }
  print(as.character(x))
  message("Alleles:", as.character(allele.names(x)), "\n")
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
#' @examples
#' \donttest{
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
