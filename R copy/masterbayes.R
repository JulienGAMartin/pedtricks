################################################################################
# genotype.list from MasterBayes

"genotype.list"<-function(G, marker.type="MSW"){
gens<-list()
if(marker.type=="MSC" | marker.type=="SNP" | marker.type=="MSW"){
  if(length(G[1,])%%2!=0){stop("Genotypes have odd number of columns")}
  for(i in 1:(length(G[1,])/2)){
    gens[[i]]<-genotype(as.matrix(G[,((i*2)-1):(i*2)]))
    names(gens)[i]<-names(G[i*2]) 
  }
}
if(marker.type=="AFLP"){
  for(i in 1:length(G[1,])){
    gens[[i]]<-genotypeD(as.matrix(G[,i]))
    names(gens)[i]<-names(G[i]) 
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

print.genotypeD <- function(x, ...) {
  if (!is.null(attr(x, "locus"))) {
    print(attr(x, "locus"))
  }
  print(as.character(x))
  cat("Alleles:", as.character(allele.names(x)), "\n")
  invisible(x)
}

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

"[.genotypeD" <- function(x, i, drop = FALSE) {
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


"[<-.genotypeD" <- function(x, i, value) {
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

"extractA"<-function(G,marker.type="MSW"){
## from Jarrod to avoid
## dependency on MasterBayes

   if(is.genotype(G[[1]])==FALSE & is.genotypeD(G[[1]])==FALSE){
     if("id"%in%colnames(G)){
       G<-G[,-which(colnames(G)=="id")]
     }
     if("categories"%in%colnames(G)){
       G<-G[,-which(colnames(G)=="categories")]
     }
     G<-genotype.list(G,marker.type=marker.type)
   }

  A<-lapply(G, function(x){summary(x)$allele.freq[,"Proportion"][1:nallele(x)]})
  if(marker.type=="AFLP" | is.genotypeD(G[[1]])){
  A<-lapply(A, function(x){
       if(any(x==1)){
         warning("some loci are monomorphic: setting allele frequencies to 0.9999 and 0.0001")
         x[which(x==1)]<-0.9999
         x[which(x==0)]<-0.0001
       }else{
         x<-x
       }
       x
     })
  }
  A
}